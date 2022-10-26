      SUBROUTINE MEDINIT(FILE,id,etam,mass,NJOB)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF 
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
C--Stores previous variables for GETTEMP (optimization)
      COMMON/PREVRUN/BPREV,XPREV,YPREV,ZPREV,TPREV,TEMPPREV
      DOUBLE PRECISION BPREV,XPREV,YPREV,ZPREV,TPREV,TEMPPREV
C--identifier of log file
	common/logfile/logfid
	integer logfid
      
      COMMON/MAXTEMPERATURE/MAXTEMP
      DOUBLE PRECISION MAXTEMP

      DATA RAU/10./
      DATA D3/0.9d0/
      DATA ZETA3/1.2d0/
C--geometry variables
      COMMON /GEOMVAR/CREATIONPOINTS,PARTICIPANTPOINTS,
     &PARTICIPANTNUMBER
      DOUBLE PRECISION PARTICIPANTPOINTS(476,2),
     &CREATIONPOINTS(2)
      DATA CREATIONPOINTS/2*0.0/
      DATA PARTICIPANTPOINTS/952*0.0/
      INTEGER PARTICIPANTNUMBER

C--local variables
      INTEGER I,LUN,POS,IOS,id,mass,TRIMLEN,NJOB
	double precision etam,getimpact,FINDMAX
      CHARACTER*100 BUFFER,LABEL
	CHARACTER*80 FILE
	character firstchar
	logical fileexist

	etamax2 = etam
	logfid = id

      IOS=0
      LUN=77
      npart_centr = 0.d0
      PARTICIPANTNUMBER = 0

C--default settings (we added SPECIES)
      TAUI=0.6d0
      TI=0.36d0
      TC=0.17d0
      WOODSSAXON=.TRUE.
      CENTRMIN=0.d0
      CENTRMAX=10.d0
      NF=3
      A=mass
      N0=0.17d0
      D=0.54d0
      SIGMANN=6.2
	MDFACTOR=0.45d0
	MDSCALEFAC=0.9d0
	boost = .true.
      SPECIES = "Pb"

C--read settings from file
	write(logfid,*)
	inquire(file=FILE,exist=fileexist)
	if(fileexist)then
        write(logfid,*)'Reading medium parameters from ',FILE
        OPEN(unit=LUN,file=FILE,status='old',err=10)
	  do 20 i=1,1000
          READ(LUN, '(A)', iostat=ios) BUFFER
	    if (ios.ne.0) goto 30
	    firstchar = buffer(1:1)
	    if (firstchar.eq.'#') goto 20
          POS=SCAN(BUFFER,' ')
          LABEL=BUFFER(1:POS)
          BUFFER=BUFFER(POS+1:)
          IF (LABEL=="TAUI")THEN
            READ(BUFFER,*,IOSTAT=IOS) TAUI
          ELSE IF (LABEL=="TI") THEN
            READ(BUFFER,*,IOSTAT=IOS) TI
          ELSE IF (LABEL=="TC") THEN
            READ(BUFFER,*,IOSTAT=IOS) TC
          ELSE IF (LABEL=="WOODSSAXON") THEN
            READ(BUFFER,*,IOSTAT=IOS) WOODSSAXON
          ELSE IF (LABEL=="CENTRMIN") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMIN
          ELSE IF (LABEL=="CENTRMAX") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMAX
          ELSE IF (LABEL=="NF") THEN
            READ(BUFFER,*,IOSTAT=IOS) NF
          ELSE IF (LABEL=="N0") THEN
            READ(BUFFER,*,IOSTAT=IOS) N0
          ELSE IF (LABEL=="D") THEN
            READ(BUFFER,*,IOSTAT=IOS) D
          ELSE IF (LABEL=="SIGMANN") THEN
            READ(BUFFER,*,IOSTAT=IOS) SIGMANN
          ELSE IF (LABEL=="MDFACTOR") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDFACTOR
          ELSE IF (LABEL=="MDSCALEFAC") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDSCALEFAC
          ELSE IF (LABEL=="SPECIES") THEN
            READ(BUFFER,*,IOSTAT=IOS) SPECIES
	    else
	      write(logfid,*)'unknown label ',label
	    endif
 20	  continue
 30	  close(LUN,status='keep')
	  write(logfid,*)'...done'
	  goto 40

 10     write(logfid,*)'Could not open medium parameter file, '//
     &	'will run with default settings.'

	else
	  write(logfid,*)'No medium parameter file found, '//
     &	'will run with default settings.'
	endif

 40   write(logfid,*)'using parameters:'
      write(logfid,*)'TAUI       =',TAUI
      write(logfid,*)'TI         =',TI
      write(logfid,*)'TC         =',TC
      write(logfid,*)'WOODSSAXON =',WOODSSAXON
      write(logfid,*)'CENTRMIN   =',CENTRMIN
      write(logfid,*)'CENTRMAX   =',CENTRMAX
      write(logfid,*)'NF         =',NF
      write(logfid,*)'A          =',A
      write(logfid,*)'N0         =',N0
      write(logfid,*)'D          =',D
      write(logfid,*)'SIGMANN    =',SIGMANN
      write(logfid,*)'MDFACTOR   =',MDFACTOR
      write(logfid,*)'MDSCALEFAC =',MDSCALEFAC
      write(logfid,*)'SPECIES    =',SPECIES
	write(logfid,*)
	write(logfid,*)
	write(logfid,*)

      SIZEC = TRIMLEN(SPECIES)
      
C-- NECESSARY CALLS TO INITIALIZE MEDIUM
      CALL SET_SEED(NJOB) !Sets the seed for the C++ function
      CALL CALCULATE_NPART_AT_CENTR(SPECIES,SIZEC,SIGMANN*10.,
     &NPART_CENTR)


      CALL CHECKFILE !Checks if the file relating b(centr) exists. Creates it if not

C-- Determination of the interval of b's for the centrality interval provided      
      bmin = getimpact(CENTRMIN)
      bmax = getimpact(CENTRMAX)

C-- Random initialization of PREVRUN parameters
      BPREV = 100.
      XPREV = 100.
      YPREV = 100.
      ZPREV = 100.
      TPREV = 100.
      TEMPPREV = 100.

C-- (PROBABLY) NOT NECESSARY
      DO 69 WHILE (PARTICIPANTNUMBER .EQ. 0)
        CALL MEDNEXTEVT
        CALL temp_glauber_medium_simple(PARTICIPANTPOINTS,
     &  PARTICIPANTNUMBER,CREATIONPOINTS,BREAL,SPECIES,SIZEC,
     &  SIGMANN*10.)
  69  CONTINUE 
      MAXTEMP = FINDMAX()

      END

C--Generation of collision via Glauber MC (added by "GAJOS DO GLAUBER")
      SUBROUTINE EVENT_GENERATION
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF 
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--geometry variables
      COMMON /GEOMVAR/CREATIONPOINTS,PARTICIPANTPOINTS,
     &PARTICIPANTNUMBER
      DOUBLE PRECISION PARTICIPANTPOINTS(476,2),
     &CREATIONPOINTS(2)
      INTEGER PARTICIPANTNUMBER
      COMMON/MAXTEMPERATURE/MAXTEMP
      DOUBLE PRECISION MAXTEMP
      DOUBLE PRECISION FINDMAX
      
      PARTICIPANTNUMBER=0
      DO 69 WHILE (PARTICIPANTNUMBER .EQ. 0)
        CALL MEDNEXTEVT
        CALL temp_glauber_medium_simple(PARTICIPANTPOINTS,
     &  PARTICIPANTNUMBER,CREATIONPOINTS,BREAL,SPECIES,SIZEC,
     &  SIGMANN*10.)
  69  CONTINUE 
      MAXTEMP = FINDMAX()
      END


C-- Checks if Centrality file exists, and creates one if it doesn't (added by "GAJOS DO GLAUBER")
      SUBROUTINE CHECKFILE
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--local variables
      integer i,j
      LOGICAL file_exists
      CHARACTER suffix*21,filename*(SIZEC+21)

      file_exists = .false.
      suffix = "--Centrality_vs_b.dat"
      DO 10 I = 1, SIZEC
        FILENAME(I:I) = SPECIES(I:I)
  10  CONTINUE
      DO 25 I=SIZEC+1,SIZEC+21
        J=I-SIZEC
        FILENAME(I:I) = SUFFIX(J:J)
  25  CONTINUE

      INQUIRE(FILE= FILENAME , EXIST=file_exists)
      IF(.NOT. file_exists) THEN
        CALL mednextevt_c(SPECIES,sizec,SIGMANN*10.,
     &number_of_lines)

        INQUIRE(FILE = FILENAME, EXIST = file_exists) 
        IF(.NOT. file_exists) THEN !file was not created
          WRITE(*,*) FILENAME // " WAS NOT CREATED!"
          STOP !program is imediately stopped
        END IF
      ELSE
        number_of_lines = 0
        OPEN(unit=98, FILE = FILENAME)
        DO
          READ(98,*,END=15)
          number_of_lines = number_of_lines + 1
        END DO
  15    CLOSE(98)
        number_of_lines = number_of_lines - 2
      ENDIF
      END


c-    Obtains the value of BREAL (added by "GAJOS DO GLAUBER")
      SUBROUTINE MEDNEXTEVT
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF  
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 

C--local variables
      DOUBLE PRECISION PYR

      breal = sqrt(PYR(0)*(bmax**2-bmin**2)+bmin**2)
      END


C--Relates a value of b with a centrality -- interpolation of file (added by "GAJOS DO GLAUBER")
      double precision function getimpact(ctemp)
      implicit none

C--medium parameters 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES

C--local variables
      integer i,j,LOGFINAL,chosen_position
      DOUBLE PRECISION b1(number_of_lines),c1(number_of_lines),
     &ctemp,btemp
      LOGICAL file_exists
      CHARACTER suffix*21,filename*(SIZEC+21)

      suffix = "--Centrality_vs_b.dat"
      DO 13 I = 1, SIZEC
        FILENAME(I:I) = SPECIES(I:I)
  13  CONTINUE
      DO 18 I=SIZEC+1,SIZEC+21
        J=I-SIZEC
        FILENAME(I:I) = SUFFIX(J:J)
  18  CONTINUE

      LOGFINAL = 99
      OPEN(UNIT=LOGFINAL,FILE= FILENAME)
C -- skip over first two lines
      read(LOGFINAL,*) 
      read(LOGFINAL,*)
C-- Stores the file points in 2 arrays
      do 20 i=1,number_of_lines
        read(LOGFINAL,*) b1(i), c1(i)
20    CONTINUE
      
      do 30 j=1,number_of_lines-1
        IF(c1(j).LE.CTEMP .AND. c1(j+1).GT.CTEMP) THEN
          chosen_position = j
          goto 11
        ENDIF
30    CONTINUE
11    btemp=0
      IF(c1(chosen_position) .NE. c1(chosen_position+1)) THEN
        btemp = (1./(c1(chosen_position)-c1(chosen_position+1)))
     &        *((b1(chosen_position)-b1(chosen_position+1))
     &        *ctemp + b1(chosen_position+1)*c1(chosen_position)
     &        - b1(chosen_position)*c1(chosen_position+1))
      ELSE 
        btemp = (0.5)*(b1(chosen_position+1)+b1(chosen_position))
      ENDIF

      getimpact=btemp
      CLOSE(LOGFINAL)
      end


      INTEGER FUNCTION TRIMLEN(STRING) !Returns length of string ignoring trailing blanks
      CHARACTER*(*) STRING
      INTEGER I
      DO 15, I = LEN(STRING), 1, -1
        IF(STRING(I:I) .NE. '') GO TO 20
  15  CONTINUE
  20  TRIMLEN = I
      END FUNCTION


C--   Chooses collision point for jet creation (modified by "GAJOS DO GLAUBER")
      SUBROUTINE PICKVTX(X,Y)
      IMPLICIT NONE
C--geometry variables
      COMMON /GEOMVAR/CREATIONPOINTS,PARTICIPANTPOINTS,
     &PARTICIPANTNUMBER
      DOUBLE PRECISION PARTICIPANTPOINTS(476,2),
     &CREATIONPOINTS(2)
      INTEGER PARTICIPANTNUMBER
C--local variables
      DOUBLE PRECISION X,Y

      X=CREATIONPOINTS(1)
      Y=CREATIONPOINTS(2)

      END


	SUBROUTINE SETB(BVAL)
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF  
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
C--local variables     
	DOUBLE PRECISION BVAL
	BREAL=BVAL
	END



      SUBROUTINE GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF  
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
C--internal medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--function calls
      DOUBLE PRECISION GETTEMP,GETMD,GETMOM,GETMS
C--identifier of log file
	common/logfile/logfid
	integer logfid
C--local variables
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E,MD,TEMP
      INTEGER TYPE
      DOUBLE PRECISION R,PYR,pmax,wt,tau,theta,phi,pi,p,ys,pz2,e2
      DATA PI/3.141592653589793d0/

      R=PYR(0)
      IF(R.LT.(2.*12.*NF*D3/3.)/(2.*12.*NF*D3/3.+3.*16.*ZETA3/2.))THEN
         TYPE=2
         MS=GETMS(X,Y,Z,T)
      ELSE
         TYPE=21
         MS=GETMD(X,Y,Z,T)
      ENDIF
      TEMP=GETTEMP(X,Y,Z,T)
	tau=sqrt(t**2-z**2)
	if (boost) then
  	  ys = 0.5*log((t+z)/(t-z))
	else
	  ys = 0.d0
	endif
	pmax = 10.*temp

      IF(TEMP.LT.1.D-2)THEN
       write(logfid,*)'asking for a scattering centre without medium:'
       write(logfid,*)'at (x,y,z,t)=',X,Y,Z,T
       write(logfid,*)'making one up to continue but '//
     &	'something is wrong!'
       TYPE=21
       PX=0.d0
       PY=0.d0
       PZ=0.d0
       MS=GETMS(0.d0,0.d0,0.d0,0.d0)
       MD=GETMD(0.d0,0.d0,0.d0,0.d0)
       E=SQRT(PX**2+PY**2+PZ**2+MS**2)
       RETURN
      ENDIF

 10	p = pyr(0)**0.3333333*pmax
	E2 = sqrt(p**2+ms**2)
	if (type.eq.2) then
	  wt = (exp(ms/temp)-1.)/(exp(E2/temp)-1.)
	else
	  wt = (exp(ms/temp)+1.)/(exp(E2/temp)+1.)
	endif
	if (wt.gt.1.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (wt.lt.0.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (pyr(0).gt.wt) goto 10
	phi = pyr(0)*2.*pi
	theta = -acos(2.*pyr(0)-1.)+pi
	px  = p*sin(theta)*cos(phi)
	py  = p*sin(theta)*sin(phi)
	pz2 = p*cos(theta)
	E   = cosh(ys)*E2 + sinh(ys)*pz2
	pz  = sinh(ys)*E2 + cosh(ys)*pz2
      END


      SUBROUTINE AVSCATCEN(X,Y,Z,T,PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
	double precision x,y,z,t,px,py,pz,e,getms,m,ys
	if (boost) then
  	  ys = 0.5*log((t+z)/(t-z))
	  if ((z.eq.0.d0).and.(t.eq.0.d0)) ys =0.d0
	  if (ys.gt.etamax2) ys=etamax2
	  if (ys.lt.-etamax2) ys=-etamax2
	else
	  ys = 0.d0
	endif
	m  = getms(x,y,z,t)
	e  = m*cosh(ys)
	px = 0.d0
	py = 0.d0
	pz = m*sinh(ys)
	end


      SUBROUTINE maxscatcen(PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
	double precision px,py,pz,e,getmsmax,m,ys
	if (boost) then
  	  ys = etamax2
	else
	  ys = 0.d0
	endif
	m  = getmsmax()
	e  = m*cosh(ys)
	px = 0.d0
	py = 0.d0
	pz = m*sinh(ys)
	end
	


      DOUBLE PRECISION FUNCTION GETMD(X1,Y1,Z1,T1)
      IMPLICIT NONE
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION X1,Y1,Z1,T1,GETTEMP
      GETMD=MDSCALEFAC*3.*GETTEMP(X1,Y1,Z1,T1)
      GETMD=MAX(GETMD,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMS(X2,Y2,Z2,T2)
      IMPLICIT NONE
      DOUBLE PRECISION X2,Y2,Z2,T2,GETMD
      GETMS=GETMD(X2,Y2,Z2,T2)/SQRT(2.)
      END



      DOUBLE PRECISION FUNCTION GETNEFF(X3,Y3,Z3,T3)
      IMPLICIT NONE
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF 
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--   local variables
      DOUBLE PRECISION X3,Y3,Z3,T3,PI,GETTEMP,tau,cosheta
      DATA PI/3.141592653589793d0/
	tau = sqrt(t3**2-z3**2)
	cosheta = t3/tau
      GETNEFF=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *GETTEMP(X3,Y3,Z3,T3)**3/PI**2
	getneff = getneff/cosheta
      END
      
      
C--Determines the temperature at a given spacetime point (modified by "GAJOS DO GLAUBER")
      DOUBLE PRECISION FUNCTION GETTEMP(X4,Y4,Z4,T4)
      IMPLICIT NONE
C--Stores previous variables for GETTEMP (optimization)
      COMMON/PREVRUN/BPREV,XPREV,YPREV,ZPREV,TPREV,TEMPPREV
      DOUBLE PRECISION BPREV,XPREV,YPREV,ZPREV,TPREV,TEMPPREV
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF 
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
      DOUBLE PRECISION X4,Y4,Z4,T4,TAU,NPART,EPS0,EPSIN,TEMPIN,PI,
     &NTHICK,ys
      LOGICAL SAMEEVENT
      DATA PI/3.141592653589793d0/

      SAMEEVENT = (X4.EQ.XPREV).AND.(Y4.EQ.YPREV).AND.(Z4.EQ.ZPREV)
     &.AND.(T4.EQ.TPREV).AND.(BREAL.EQ.BPREV) 
      IF(SAMEEVENT) then
        GETTEMP = TEMPPREV
c      ENDIF
      ELSE

        GETTEMP=0.D0
  
        IF(ABS(Z4).GT.T4)RETURN
  
        TAU=SQRT(T4**2-Z4**2)
C--check for overlap region --- MODIFIED (DEPENDENCE OF NTHICK MUST BE ELIMINATED)
        IF(NPART(X4,Y4) .EQ. 0.d0) RETURN
  
	  ys = 0.5*log((t4+z4)/(t4-z4))
	  if (abs(ys).gt.etamax2) return
C--determine initial temperature at transverse position
        IF(WOODSSAXON)THEN
           EPS0=(16.*8.+7.*2.*6.*NF)*PI**2*TI**4/240.
           EPSIN=EPS0*NPART(X4,Y4)
           TEMPIN=((EPSIN*240./(PI**2*(16.*8.+7.*2.*6.*NF)))
     &        **0.25)/NPART_CENTR
        ELSE
           TEMPIN=TI
        ENDIF
C--calculate temperature if before initial time
        IF(TAU.LE.TAUI)THEN
	   GETTEMP=TEMPIN*TAU/TAUI
        ELSE
C--evolve temperature
         GETTEMP=TEMPIN*(TAUI/TAU)**0.3333
        ENDIF
        IF(GETTEMP.LT.TC) GETTEMP=0.d0
        BPREV = BREAL
        XPREV = X4
        YPREV = Y4
        ZPREV = Z4
        TPREV = T4
        TEMPPREV = GETTEMP
      ENDIF
 
      END

C--Set Search region to look for maximum
      SUBROUTINE SET_SEARCH_REGION(XMIN,XMAX,YMIN,YMAX)
      IMPLICIT NONE 
      COMMON /GEOMVAR/CREATIONPOINTS,PARTICIPANTPOINTS,
     &PARTICIPANTNUMBER
      DOUBLE PRECISION PARTICIPANTPOINTS(476,2),
     &CREATIONPOINTS(2)
      INTEGER PARTICIPANTNUMBER

C-- Local variables
      DOUBLE PRECISION XMIN,XMAX,YMIN,YMAX,EPSILON
      INTEGER i
      XMIN = PARTICIPANTPOINTS(1,1)
      YMIN = PARTICIPANTPOINTS(1,2)
      XMAX = PARTICIPANTPOINTS(1,1)
      YMAX = PARTICIPANTPOINTS(1,2)
      EPSILON = 0.4
      

      DO 96 I=2,PARTICIPANTNUMBER
        IF(PARTICIPANTPOINTS(I,1).LT.XMIN) THEN
          XMIN = PARTICIPANTPOINTS(I,1)
        ELSE IF(PARTICIPANTPOINTS(I,1).GT.XMAX) THEN
          XMAX = PARTICIPANTPOINTS(I,1)
        ENDIF

        IF(PARTICIPANTPOINTS(I,2).LT.YMIN) THEN
          YMIN = PARTICIPANTPOINTS(I,2)
        ELSE IF(PARTICIPANTPOINTS(I,2).GT.YMAX) THEN
          YMAX = PARTICIPANTPOINTS(I,2)
        ENDIF
  96  CONTINUE
      
      XMIN = XMIN - EPSILON
      YMIN = YMIN - EPSILON
      XMAX = XMAX + EPSILON
      YMAX = YMAX + EPSILON

      END

C-- Finds the maximum
      DOUBLE PRECISION FUNCTION FINDMAX()
      IMPLICIT NONE
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES


C-- Local Variables
      INTEGER XPOINTS,YPOINTS,IX,IY
      DOUBLE PRECISION XSTEP,YSTEP,GETTEMP,TEMPRESULT
      DOUBLE PRECISION XMIN,XMAX,YMIN,YMAX,DBLE
      DOUBLE PRECISION XTEMP,YTEMP,XFINAL,YFINAL

      FINDMAX = 0.
      XFINAL = 0.
      YFINAL = 0.

      XPOINTS = 200
      YPOINTS = 200
      
      CALL SET_SEARCH_REGION(XMIN,XMAX,YMIN,YMAX)

      
      XSTEP = (XMAX-XMIN)/DBLE(XPOINTS-1)
      YSTEP = (YMAX-YMIN)/DBLE(YPOINTS-1)

      DO 95 IX=0,XPOINTS-1
        XTEMP = XMIN + IX*XSTEP
        DO 94 IY=0,YPOINTS-1
          YTEMP = YMIN + IY*YSTEP
          TEMPRESULT=GETTEMP(XTEMP,YTEMP,0.D0,TAUI)
          IF(TEMPRESULT.GT.FINDMAX) THEN
            FINDMAX = TEMPRESULT
            XFINAL = XTEMP
            YFINAL = YTEMP
          ENDIF
  94    CONTINUE
  95  CONTINUE      
      END


      DOUBLE PRECISION FUNCTION GETTEMPMAX()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF 
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
      COMMON/MAXTEMPERATURE/MAXTEMP
      DOUBLE PRECISION MAXTEMP
      double precision gettemp

C--function call
      GETTEMPMAX = MAXTEMP
c      GETTEMPMAX = GETTEMP(0.D0,0.D0,0.D0,TAUI)
      
      END



      DOUBLE PRECISION FUNCTION GETMDMAX()
      IMPLICIT NONE
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
      GETMDMAX=MDSCALEFAC*3.*GETTEMPMAX()
      GETMDMAX=MAX(GETMDMAX,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMDMIN()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
	GETMDMIN=MDSCALEFAC*3.*TC
      GETMDMIN=MAX(GETMDMIN,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMSMAX()
      IMPLICIT NONE
      DOUBLE PRECISION GETMDMAX,SQRT
      GETMSMAX=GETMDMAX()/SQRT(2.D0)
      END



	DOUBLE PRECISION FUNCTION GETNATMDMIN()
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF 
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC,PI
      DATA PI/3.141592653589793d0/
C--local variables
	DOUBLE PRECISION T,GETMDMIN
	T=GETMDMIN()/(MDSCALEFAC*3.)
      GETNATMDMIN=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *T**3/PI**2
	END



	DOUBLE PRECISION FUNCTION GETLTIMEMAX()
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF  
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--function call
      DOUBLE PRECISION GETTEMPMAX
	GETLTIMEMAX=TAUI*(GETTEMPMAX()/TC)**3*cosh(etamax2)
	END



      DOUBLE PRECISION FUNCTION GETNEFFMAX()
      IMPLICIT NONE
      COMMON/MEDPARAM/BMAX,BMIN,CENTRMIN,CENTRMAX,BREAL,CENTR,
     &RAU,NF 
      INTEGER NF
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     &BMAX,BMIN 
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,SPECIES,SIZEC,
     &number_of_lines,npart_centr
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN,npart_centr
      INTEGER A,SIZEC,number_of_lines
      LOGICAL WOODSSAXON
      CHARACTER*80 SPECIES
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--   local variables
      DOUBLE PRECISION PI,GETTEMPMAX
      DATA PI/3.141592653589793d0/
      GETNEFFMAX=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *GETTEMPMAX()**3/PI**2
      END
      
      
C--   Density of participants (modified by "GAJOS DO GLAUBER")
      DOUBLE PRECISION FUNCTION NPART(XX,YY)
      IMPLICIT NONE
C--geometry variables
      COMMON /GEOMVAR/CREATIONPOINTS,PARTICIPANTPOINTS,
     &PARTICIPANTNUMBER
      DOUBLE PRECISION PARTICIPANTPOINTS(476,2),
     &CREATIONPOINTS(2)
      INTEGER PARTICIPANTNUMBER
C--local variables
      DOUBLE PRECISION XX,YY,SIGMA,GAUSSIAN2D
      INTEGER I,limit

      SIGMA = 0.4 !standard deviation of gaussian
      NPART = 0.
C--Only participant points that are limit*sigma distance away from (XX,YY) are considered
C--(if limit = 0, then all participant points are considered)
      limit = 4 
      IF(PARTICIPANTNUMBER.EQ.0) RETURN

      DO 33 I=1,PARTICIPANTNUMBER

C     Equal standard deviation (sigmax=sigmay=sigma)
      IF((((xx-PARTICIPANTPOINTS(I,1))**2 + 
     &(yy-PARTICIPANTPOINTS(I,2))**2).LE.
     &(limit*limit*SIGMA*SIGMA)).OR.limit.EQ.0) Then
        NPART = NPART + GAUSSIAN2D(XX,YY,PARTICIPANTPOINTS(I,1),
     &  PARTICIPANTPOINTS(I,2),SIGMA)
      ENDIF

   33 CONTINUE
      
      END

C--   2D Gaussian formula, normalized to unit (added for NPART function by "GAJOS DO GLAUBER")
      DOUBLE PRECISION FUNCTION GAUSSIAN2D(X,Y,XC,YC,SIGMA)
      IMPLICIT NONE
C--local variables
      DOUBLE PRECISION X,Y,XC,YC,SIGMA,PI
      DATA PI/3.141592653589793d0/

      GAUSSIAN2D = (1./(2.*PI*SIGMA*SIGMA))*
     &EXP(-(((X-XC)*(X-XC)+(Y-YC)*(Y-YC))/(2.*SIGMA*SIGMA)))

      END