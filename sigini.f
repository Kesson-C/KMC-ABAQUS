c                    
c                    user subroutine sigini

      SUBROUTINE SIGINI(SIGMA,COORDS,NTENS,NCRDS,NOEL,NPT,LAYER,
     1                  KSPT,LREBAR,NAMES)
c-----------------------------------------------------------------------------------       
c     This subroutine is used to define the initial stresses which are 
c     calculated from kinetic Monte Carlo simulation method (KMC) for 
c     all elements in the model. 
c
c     This subrouitine will be called for user-subroutine-defined initial 
c     stress fields at particular material points (these are the effective 
c     stress values for soils analysis). In addition, this subroutine will
c     be called at each iteration point in elements.
c
c     Input variables  
c
c       COORDS: An array containing the initial coordinates of this point.
c       NTENS: Number of stresses to be defined, which depends on the element type.
c       NCRDS: Number of coordinates.
c       NOEL: Element number.
c       NPT: Integration point number in the element.
c       LAYER: Layer number (for composite shells and layered solids).
c       KSPT: Section point number within the current layer.
c
c     Working variables:
c       NTELE: number of total elements in the model
c       NRELE: number of elements with residual stresses in the model = NTELE
c-----------------------------------------------------------------------------------
c
      include 'aba_param.inc'
c                    
      character*80 cmname
c                    
      dimension SIGMA(NTENS),COORDS(NCRDS)
      integer ntele,nrthe
      real,save :: RES(800000,3)
      real Elas,nu,DP
      integer, save :: nsig
c     real, save, allocatable :: RES(:,:)
c                    
c     (1) Open the input file which included the eigrnstrains
c
      open (unit=31,file="C:\Temp\Res_strain.txt")
c      
c     (2) Read residual strains from KMC results and calculate
c         the corresponding residual stresses in order of S11, S12, S22
c         since the residual strains are given in order of E11, E12, E33
c         ( dynamic array to store the stress information )
c
      read (31,*) Elas,nu
      read (31,*) ntele
c      allocate ( RES(ntele,3) )
      if (nsig == 0 ) then
          do ie = 1, ntele
             read (31,*) k,(RES(ie,j),j=1,3)
          enddo   
          DP = Elas/(1.0+nu)/(1.0-2*nu)
          do ie = 1, ntele
             RES(ie,1) = DP*( (1-nu)*RES(ie,1)+nu*RES(ie,3) ) 
             RES(ie,2) = DP*(1-2*nu)*RES(ie,2)
             RES(ie,3) = DP*( (1-nu)*RES(ie,3)+nu*RES(ie,1) )
          enddo                                      
      endif       
c
c     (3) Assign residual stress to each element 
c         Note : The residual stress vector in ABAQUS resuires 
c                in order with S11, S22, S12 
c      
      do ie = 1,ntele
         if (ie .eq. noel) then    
             SIGMA(1) = RES(NOEL,1)
             SIGMA(2) = RES(NOEL,3)
             SIGMA(3) = RES(NOEL,2)    
             goto 20
         endif    
      enddo
c 
 20   close (unit=31)
c      close (unit=32)
      nsig = nsig + 1
c            
      return
      end

