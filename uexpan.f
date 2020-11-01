c                    
c                    user subroutine uexpan

      subroutine uexpan(expan,dexpandt,temp,time,dtime,predef,dpred,
     $     statev,cmname,nstatv,noel)
c       
c     This subroutine is used to define incremental thermal strains 
c     (eigenstrain) as functions of temperature, predefined field 
c     variables, and state variables.
c     This subrouitine is called at all integration points of elements
c     for which the material behavior definition contains user-subroutine-
c     defined thermal expansion.
c
      include 'aba_param.inc'
c                    
      character*80 cmname
c                    
      dimension expan(*),dexpandt(*),temp(2),time(2),predef(*),
     $     dpred(*),statev(nstatv)
      real strain
      real,save :: Str(800000,3)      
      integer NTele,ie,nthe
c                    
c     (1) Open the input file which included the eigrnstrains
c
      open (unit=21,file="C:\Temp\Eigenstrain.txt")
c     open (unit=22,file="C:\Temp\PY_Example1_001.chk")
c     write (22,*) "Element ID from ABAQUS: ",noel
c      
c     (2) Read input file and assign eigenstrain to the corresponding element
c
      if (noel == 1 ) then
          read (21,*) ntele
          do ie = 1, ntele
             read (21,*) k,(Str(ie,j),j=1,3)
          enddo   
      endif    
c
      expan(1) = Str(noel,1)
      expan(2) = Str(noel,3)      
      expan(3) = Str(noel,2)
c     write (22,*) noel,expan(3)
c 
      close (unit=21)      
c            
      return
      end

