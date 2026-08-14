"""Optional: view the wheel in the Autodesk Platform Services Viewer, inside Colab.

The built-in Plotly viewer draws the deck's own triangles and is verified vertex-for-
vertex against the .inp. This is an alternative that hands the geometry to Autodesk's
renderer instead -- better shading, section planes, model tree -- at the cost of a
round trip through their cloud.

What a view costs, every time
-----------------------------
1. an OAuth token from your APS app (client id + secret)
2. the STEP/STL uploaded to an OSS bucket
3. a Model Derivative job translating it to SVF2  <- this is what consumes flex tokens
4. polling until the translation finishes         <- minutes on a large model
5. the viewer embedded in the output cell with a live token

So it is not a substitute for the local preview during iteration. Use it when you want
a presentation-quality view of a model you have already settled on.

Honesty about this module
-------------------------
The endpoint paths and payloads here are written from knowledge of the APS REST API,
**not** verified against the live service -- the documentation is JavaScript-rendered
and could not be read, and there were no credentials to test with. Every call therefore
checks its status code and prints the server's own error body, so the first real run
tells you exactly which call is wrong rather than failing opaquely.

Nothing here is imported by the rest of semgrit. If you never call it, it costs nothing.
"""

from __future__ import annotations

import base64
import json
import os
import time
from dataclasses import dataclass, field
from typing import Optional

BASE = "https://developer.api.autodesk.com"
AUTH = BASE + "/authentication/v2/token"
OSS = BASE + "/oss/v2/buckets"
MD = BASE + "/modelderivative/v2/designdata"

# A translation is billed per job and a big multi-body STEP is the expensive case, so
# refuse rather than surprise the user with a bill. Both are overridable.
DEFAULT_MAX_MB = 100.0
DEFAULT_MAX_BODIES = 2000


class APSError(RuntimeError):
    pass


@dataclass
class APSConfig:
    client_id: str = ""
    client_secret: str = ""
    bucket_key: str = ""
    """Must be globally unique and lowercase. Blank = derived from the client id."""
    region: str = "US"                 # 'US' | 'EMEA'
    max_upload_mb: float = DEFAULT_MAX_MB
    max_bodies: int = DEFAULT_MAX_BODIES
    poll_seconds: float = 5.0
    timeout_minutes: float = 20.0

    def resolved_bucket(self) -> str:
        if self.bucket_key:
            return self.bucket_key.lower()
        # Bucket keys share one global namespace, so tie it to the app that owns it.
        return ("semgrit-" + self.client_id.lower())[:120]


def _req(method, url, *, headers=None, data=None, json_body=None, expect=(200, 201)):
    """One HTTP call, with the server's error body surfaced rather than swallowed."""
    import requests

    r = requests.request(method, url, headers=headers or {}, data=data,
                         json=json_body, timeout=300)
    if r.status_code not in expect:
        body = r.text[:600]
        raise APSError("%s %s -> HTTP %d\n%s" % (method, url, r.status_code, body))
    if r.content and r.headers.get("content-type", "").startswith("application/json"):
        return r.json()
    return r


def get_token(cfg: APSConfig, scopes: str = "data:read data:write data:create "
                                            "bucket:create bucket:read viewables:read"):
    """Two-legged OAuth. Returns (access_token, expires_in_seconds)."""
    if not cfg.client_id or not cfg.client_secret:
        raise APSError("no APS client id/secret. Create an app at aps.autodesk.com "
                       "(Custom Integration, with Model Derivative and Data "
                       "Management enabled) and paste its credentials -- an Autodesk "
                       "account email and password will not work here.")
    basic = base64.b64encode(
        ("%s:%s" % (cfg.client_id, cfg.client_secret)).encode()).decode()
    out = _req("POST", AUTH,
               headers={"Authorization": "Basic " + basic,
                        "Content-Type": "application/x-www-form-urlencoded"},
               data={"grant_type": "client_credentials", "scope": scopes})
    return out["access_token"], out.get("expires_in", 3600)


def ensure_bucket(cfg: APSConfig, token: str) -> str:
    bucket = cfg.resolved_bucket()
    h = {"Authorization": "Bearer " + token, "Content-Type": "application/json"}
    try:
        _req("GET", "%s/%s/details" % (OSS, bucket), headers=h)
        return bucket
    except APSError:
        pass
    _req("POST", OSS, headers={**h, "x-ads-region": cfg.region},
         json_body={"bucketKey": bucket, "policyKey": "transient"},
         expect=(200, 201, 409))     # 409 = someone already owns that key
    return bucket


def check_size(path: str, cfg: APSConfig) -> dict:
    """Refuse before uploading anything, not after."""
    mb = os.path.getsize(path) / 1e6
    info = {"path": path, "mb": mb, "bodies": None}
    if path.lower().endswith(".step") or path.lower().endswith(".stp"):
        # Count solid bodies cheaply: a multi-body STEP of a dressed wheel can carry
        # thousands, and Model Derivative charges for the pain. Count every B-rep
        # flavour, not just MANIFOLD_SOLID_BREP -- semgrit writes FACETED_BREP, so a
        # counter that only knew the manifold form silently reported zero bodies and
        # the cap never fired. Scan in chunks: a STEP DATA section can be one line.
        BREPS = ("FACETED_BREP", "MANIFOLD_SOLID_BREP", "BREP_WITH_VOIDS",
                 "SHELL_BASED_SURFACE_MODEL", "MANIFOLD_SURFACE_SHAPE_REPRESENTATION")
        n, tail = 0, ""
        with open(path, "r", encoding="ascii", errors="ignore") as fh:
            while True:
                chunk = fh.read(1 << 20)
                if not chunk:
                    break
                buf = tail + chunk
                for kw in BREPS:
                    n += buf.count(kw)
                tail = buf[-64:]        # keep an overlap so a split keyword is caught
        info["bodies"] = n
        if n > cfg.max_bodies:
            raise APSError(
                "%s holds %d solid bodies, over the %d cap. Translating a multi-body "
                "STEP of that size is slow and token-expensive -- send the STL "
                "instead, or cap the grit count, or raise max_bodies deliberately."
                % (os.path.basename(path), n, cfg.max_bodies))
    if mb > cfg.max_upload_mb:
        raise APSError(
            "%s is %.1f MB, over the %.1f MB cap. Upload and translation are billed, "
            "so this is refused rather than charged. Send the STL rather than the "
            "STEP, lower the grit count, or raise max_upload_mb deliberately."
            % (os.path.basename(path), mb, cfg.max_upload_mb))
    return info


def upload(cfg: APSConfig, token: str, bucket: str, path: str) -> str:
    """Signed-S3 upload. Returns the base64url URN the viewer needs."""
    import requests

    obj = os.path.basename(path)
    h = {"Authorization": "Bearer " + token}
    signed = _req("GET", "%s/%s/objects/%s/signeds3upload?minutesExpiration=30"
                  % (OSS, bucket, obj), headers=h)
    urls = signed["urls"]
    with open(path, "rb") as fh:
        blob = fh.read()
    # One part is enough below 100 MB; the cap above keeps us there.
    put = requests.put(urls[0], data=blob, timeout=900)
    if put.status_code not in (200, 201):
        raise APSError("S3 PUT -> HTTP %d\n%s" % (put.status_code, put.text[:400]))
    done = _req("POST", "%s/%s/objects/%s/signeds3upload" % (OSS, bucket, obj),
                headers={**h, "Content-Type": "application/json"},
                json_body={"uploadKey": signed["uploadKey"]})
    object_id = done["objectId"]
    return base64.urlsafe_b64encode(object_id.encode()).decode().rstrip("=")


def translate(cfg: APSConfig, token: str, urn: str) -> None:
    _req("POST", MD + "/job",
         headers={"Authorization": "Bearer " + token,
                  "Content-Type": "application/json", "x-ads-force": "true"},
         json_body={"input": {"urn": urn},
                    "output": {"formats": [{"type": "svf2", "views": ["3d"]}]}},
         expect=(200, 201, 202))


def wait(cfg: APSConfig, token: str, urn: str, log=print) -> dict:
    h = {"Authorization": "Bearer " + token}
    deadline = time.time() + cfg.timeout_minutes * 60
    last = None
    while time.time() < deadline:
        m = _req("GET", "%s/%s/manifest" % (MD, urn), headers=h)
        status, prog = m.get("status"), m.get("progress", "")
        if (status, prog) != last:
            log("  translation: %s %s" % (status, prog))
            last = (status, prog)
        if status == "success":
            return m
        if status in ("failed", "timeout"):
            raise APSError("translation %s:\n%s" % (status, json.dumps(m)[:800]))
        time.sleep(cfg.poll_seconds)
    raise APSError("translation still running after %g minutes; re-run wait() with the "
                   "same URN rather than translating again" % cfg.timeout_minutes)


def viewer_html(urn: str, token: str, height: int = 640) -> str:
    """HTML for a Colab output cell. The token is short-lived and scoped read-only."""
    ver = "7.*"
    return """
<link rel="stylesheet" href="https://developer.api.autodesk.com/modelderivative/v2/viewers/%(v)s/style.min.css" type="text/css">
<script src="https://developer.api.autodesk.com/modelderivative/v2/viewers/%(v)s/viewer3D.min.js"></script>
<div id="apsViewer" style="position:relative;width:100%%;height:%(h)dpx;"></div>
<script>
  var opts = {
    env: 'AutodeskProduction2',
    api: 'streamingV2',
    getAccessToken: function (cb) { cb('%(t)s', 3600); }
  };
  Autodesk.Viewing.Initializer(opts, function () {
    var div = document.getElementById('apsViewer');
    var v = new Autodesk.Viewing.GuiViewer3D(div);
    v.start();
    Autodesk.Viewing.Document.load('urn:%(u)s', function (doc) {
      v.loadDocumentNode(doc, doc.getRoot().getDefaultGeometry());
    }, function (code, msg) {
      div.innerHTML = '<pre>viewer could not load the model: ' + code + ' ' + msg +
                      '</pre>';
    });
  });
</script>
""" % {"v": ver, "h": height, "t": token, "u": urn}


def probe_html(height: int = 120) -> str:
    """Does Colab's sandboxed output iframe even allow the viewer library to load?

    Worth answering before creating an account: if the sandbox blocks the script or
    its XHRs, no amount of correct API code will help.
    """
    return """
<div id="apsProbe" style="font-family:monospace;font-size:13px">loading...</div>
<script src="https://developer.api.autodesk.com/modelderivative/v2/viewers/7.*/viewer3D.min.js"
        onload="document.getElementById('apsProbe').innerHTML =
                 'viewer3D.js LOADED, Autodesk.Viewing is ' +
                 (window.Autodesk && window.Autodesk.Viewing ? 'available'
                                                             : 'MISSING');"
        onerror="document.getElementById('apsProbe').innerHTML =
                 'BLOCKED: the sandbox refused to load viewer3D.js';"></script>
""" % {}


def publish(path: str, cfg: APSConfig, log=print) -> dict:
    """Size-check, upload, translate, wait. Returns {urn, token, manifest, size}."""
    size = check_size(path, cfg)
    log("file       : %s, %.1f MB%s"
        % (os.path.basename(path), size["mb"],
           ", %d solid bodies" % size["bodies"] if size["bodies"] else ""))
    token, _ = get_token(cfg)
    log("auth       : token acquired")
    bucket = ensure_bucket(cfg, token)
    log("bucket     : %s" % bucket)
    urn = upload(cfg, token, bucket, path)
    log("uploaded   : urn %s..." % urn[:28])
    translate(cfg, token, urn)
    manifest = wait(cfg, token, urn, log=log)
    log("translated : ready to view")
    return {"urn": urn, "token": token, "manifest": manifest, "size": size}
