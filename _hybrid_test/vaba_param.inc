c
c     Stand-in for Abaqus' own vaba_param.inc, so the VUMATs can be
c     compiled and exercised outside Abaqus. Abaqus/Explicit ships the
c     double precision variant when the job is submitted with
c     double=both, which is what this project's decks require.
c
      implicit real*8(a-h,o-z)
      parameter (nprecd=2)
