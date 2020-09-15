c
c     Code to replicate the results from the purely-lytic vs temperate
c     models in Knowles et al "Temperate infection in a virus–host
c     system previously known for virulent dynamics", Nature
c     Communications 11, 4626 (2020).
c
c     The code below integrates numerically Eqs.(1)-(3) in the paper,
c     and associated equations (Eqs.(4)-(7)), which can be used to model
c     the dynamics of a host i) in the absence of viruses; ii) in the
c     presence of a purely lytic virus; or iii) in the presence of a
c     temperate virus, as described in the manuscript.
c      
c     For case i) (host only), set parameter MOI0 (initial ratio
c     virus:host density) to zero. For cases with viruses, set MOI0 to
c     10. For case ii) (purely lytic virus), set ipurel to 1, and for
c     case iii) (temperate virus), set ipurel to 0.
c
c     To replicate the standard infection experiments ("ambient soup"),
c     set ipre to 0 and change the initial density of hosts Ehux0 as
c     needed (see units below).
c
c     To replicate the pre-infection treatments, set ipre to 1. The code
c     the will first generate a "natural initial condition" as in the
c     experimental setup, by letting a system with Ehux0 initial host
c     density and EhV0=MOI0*Ehux0 intial density interact for tpre
c     days. Following the experiments, the resulting system is then
c     "diluted" to achieve the desired initial density for hosts, and
c     all viruses are removed. To replicate the experiments with an
c     additional initial inoculum of viruses set ispk to 1, which adds
c     1E6 viruses per ml to the diluted system.
c      

      Program replicate_Ehux_EhV_infection_in_batch_environments
      Implicit none
      Double precision Ehux_free,EhV,EhuxI,Ehux,Ehux_deaths,t,dt,mu,
     & mu_max,mort,k,B,L,mortV,MOI,tmax,MOI0,Ehux0,EhV0,dEhuxFdt,dEhVdt,
     &dEhuxIdt,Ehux_lysis,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,N0,xaux,Vbatch,
     &EhuxF_deaths,EhuxI_deaths,tpre,muI,Ehux_OLD,mu_min,tmumin,Kcarr,
     &slope,intercept,mu_eff,tphoto,
c     Temporal observables:
     &mu_t,Ehux_t,EhuxF_t,EhV_t,EhuxI_t,points,MOI_t,Ehux_deaths_t,time,
     &net_mut,
c     Fake integer
     &it_fake
      Integer i,j,krun,Kmax,it,it_old,itmax,steps,krun_eff,ipre,ispk,
     &ilytic,itemp,ipurel
      Character*5 ck
      PARAMETER (tmax=5d1,dt=1d-4,steps=300,ispk=0,ipre=0,ipurel=0
     &,itmax=tmax/dt/steps,MOI0=1d1,tpre=1d0+2d0/24d0,Kmax=1)      
      Dimension points(itmax),mu_t(itmax),Ehux_t(itmax),EhuxF_t(itmax),
     &EhV_t(itmax),EhuxI_t(itmax),Ehux_deaths_t(itmax),time(itmax),
     &MOI_t(itmax),net_mut(itmax)

      
      N0=0.2581*900d-6          !mol/l of nitrogen in f/2 media, Ben
      Vbatch=1d0                !L
c     
c     Ehux parametrization
c
      mu_max=1.01d0             !Max. growth rate, d^-1
      mort=0d0                  !Natural mortality rate, d^-1
c
c     EhV parametrization
c
      k=10*1d-9*1d-3*60*24d0    !10*1d-9 ml/min in l/d
      B=100                     !virions per infection
      L=2d0                     !days
      mortV=1d0/3d0             !viral infectivity decline rate, d^-1

c
c     Population initialization
c      
      Ehux0=1d1*1d3             !cell/l
      EhV0=MOI0*Ehux0           !viruses/l

      Kcarr=6.5d6*1d3           !Max population density in the f/2
                                !experiments

c
c     Calculation of the parameters involved in Eq.(5); for the exact
c     densities from the manuscript, use the commented parts, but the
c     non-commented general expression enables extrapolation to any
c     initial density:
c      
      mu_min=0.75d0
      tmumin=0d0
      IF (Ehux0.eq.1d1*1d3) THEN
         tmumin=19d0
         mu_min=0.7d0
      END IF
c      IF (Ehux0.eq.1d2*1d3) tmumin=17d0
c      IF (Ehux0.eq.1d3*1d3) tmumin=14d0
c      IF (Ehux0.eq.1d4*1d3) tmumin=9d0
c      IF (Ehux0.eq.1d5*1d3) tmumin=5d0
c      IF (Ehux0.eq.1d6*1d3) tmumin=2d0
      IF (Ehux0.gt.1d1*1d3) tmumin=-1.75d0*dlog(Ehux0*1d-3)+25d0

      slope=(mu_min-mu_max)/(tmumin-2d0)
      intercept=mu_max-slope*2d0

c
c     Calculation of the time at which photosynthetic efficiency drops
c     for the uninfected case, f/2 experiments, using the general
c     expression described in Supplementary Material:
c      

      tphoto=1d0
      tphoto=-1.3d0*dlog(Ehux0*1d-3)+16d0
      IF (Ehux0.gt.3.5d7) tphoto=1d0

c     The arrays below store the dynamics/values for the variables for
c     all replicates, in case we do want to implement several replicates
c     (Kmax). Using several replicates is important if accounting for
c     small variability across replicates (e.g. by including
c     stochasticity in either parameter or dynamics).
      
      DO i=1,itmax
         time(i)=0d0
         points(i)=0d0
         mu_t(i)=0d0
         Ehux_t(i)=0d0
         EhuxF_t(i)=0d0
         EhV_t(i)=0d0
         EhuxI_t(i)=0d0
         Ehux_deaths_t(i)=0d0
         MOI_t(i)=0d0
         net_mut(i)=0d0
      END DO

      krun_eff=0
      DO krun=1,Kmax            
         
         write(ck,fmt='(i5.5)') krun

         Ehux_free=Ehux0
         EhuxI=Ehux0-Ehux_free
         EhV=EhV0

         MOI=MOI0

c
c     Preinfection:
c
         t=1d0
         itemp=1
         IF (ipre.eq.1) THEN

            Ehux_free=1d6*1d3   !cell per l
            EhV=1d1*Ehux_free   !10:1 virus to host ratio
            DO WHILE (t.le.tpre)

c
c     Lytic-temperate switch. While virus is temperate, itemp=1. When
c     stress levels trigger lysis (time represented by tphoto, see
c     Eq.(7) and explanation below in main text), ilytic turns 1. For
c     purely lytic viruses, ilytic=1 throughout the whole simulation.
c     Although, as in our experiments, this switch never happens in the
c     pre-infection component of the experiment, for coherence we used
c     the exact same code for these dynamics as for the rest of the
c     experiment below:
c            
               itemp=(1-ipurel)
               IF (t.ge.tphoto) itemp=0
               ilytic=1-itemp

c
c     Update of host growth rate:
c     
               
               mu_eff=mu_max
               IF ((t-1d0.ge.2d0).and.(mu_max.gt.0.75d0)) THEN
                  mu_eff=slope*(t-1d0)+intercept
                  IF (t-1d0.ge.tmumin) mu_eff=mu_min
               END IF
               
               mu=mu_eff*(1d0-(Ehux_free+EhuxI*itemp)/Kcarr)
               muI=mu*itemp

               Ehux_lysis=ilytic*(EhuxI/L)*dt

               EhuxI_deaths=mort*EhuxI*dt+Ehux_lysis
               EhuxF_deaths=mort*Ehux_free*dt
               Ehux_deaths=EhuxF_deaths+EhuxI_deaths
                      
               dEhuxFdt=(mu-k*EhV/Vbatch)*Ehux_free+muI*EhuxI
               dEhuxIdt=k*Ehux_free/Vbatch*EhV
               dEhVdt=-(k*Ehux/Vbatch+mortV)*EhV
            
               Ehux_free=Ehux_free+dEhuxFdt*dt-EhuxF_deaths

               EhuxI=EhuxI+dEhuxIdt*dt-EhuxI_deaths
              
               EhV=EhV+B*Ehux_lysis+dEhVdt*dt
            
               Ehux=Ehux_free+EhuxI

               MOI=EhV/Ehux

               t=t+dt
               
            END DO

            xaux=EhuxI/Ehux
            EhuxI=Ehux0*xaux

            Ehux_free=Ehux0-EhuxI
            EhV=ispk*1d6*1d3         !Spike treatment, added 1d6 viruses/ml

            Ehux=Ehux_free+EhuxI
            
         END IF

         t=1d0
         it=1
         it_fake=it
         DO WHILE (t.le.tmax)

            it_old=it          
            IF (mod(it_fake,steps+0d0).eq.0) it=it+1

            Ehux_OLD=Ehux

c
c     Lytic-temperate switch. While virus is temperate, itemp=1. When
c     stress levels trigger lysis (time represented by tphoto, see
c     Eq.(7) and explanation below in main text), ilytic turns 1. For
c     purely lytic viruses, ilytic=1 throughout the whole simulation:
c            
            itemp=(1-ipurel)
            IF (t.ge.tphoto) itemp=0
            ilytic=1-itemp

c
c     Update of host growth rate:
c            
            
            mu_eff=mu_max
            IF ((t-1d0.ge.2d0).and.(mu_max.gt.0.75d0)) THEN
               mu_eff=slope*(t-1d0)+intercept
               IF (t-1d0.ge.tmumin) mu_eff=mu_min
            END IF
               
            mu=mu_eff*(1d0-(Ehux_free+EhuxI*itemp)/Kcarr)
            muI=mu*itemp
            
            Ehux_lysis=ilytic*(EhuxI/L)*dt
            
            EhuxI_deaths=mort*EhuxI*dt+Ehux_lysis
            EhuxF_deaths=mort*Ehux_free*dt
            Ehux_deaths=EhuxF_deaths+EhuxI_deaths
                      
            dEhuxFdt=(mu-k*EhV/Vbatch)*Ehux_free+muI*EhuxI
            
            dEhuxIdt=k*Ehux_free/Vbatch*EhV
            
            dEhVdt=-(k*Ehux/Vbatch+mortV)*EhV
            
            Ehux_free=Ehux_free+dEhuxFdt*dt-EhuxF_deaths

            EhuxI=EhuxI+dEhuxIdt*dt-EhuxI_deaths
              
            EhV=EhV+B*Ehux_lysis+dEhVdt*dt
            
            Ehux=Ehux_free+EhuxI

            MOI=EhV/Ehux

            
c*******************************************
c     Calculation of observables:
c*******************************************

cx            IF (it.ne.it_old) THEN   !In case we measure only at specific times
               mu_t(it)=mu_t(it)+mu
               
               EhuxF_t(it)=EhuxF_t(it)+Ehux_free
               EhuxI_t(it)=EhuxI_t(it)+EhuxI
               Ehux_t(it)=Ehux_t(it)+Ehux
               EhV_t(it)=EhV_t(it)+EhV
               Ehux_deaths_t(it)=Ehux_deaths_t(it)+Ehux_deaths

               MOI_t(it)=MOI_t(it)+MOI

               xaux=(Ehux-Ehux_OLD)/Ehux/dt
               IF (Ehux_OLD.eq.0d0) xaux=mu !For the very first point
               net_mut(it)=net_mut(it)+xaux

               time(it)=time(it)+t
               points(it)=points(it)+1d0
cx            END IF


            t=t+dt              !d
            it_fake=anint(it_fake+1d0)

c
c     Condition for extinction, to avoid artificial (i.e. numerically
c     driven) population recovery from unrealistic low densities:
c            
            IF (Ehux.lt.1d0) THEN
               
               write(*,*) 'The total host population collapses!!',t
               GO TO 10
               
            END IF

         END DO

 10      CONTINUE

         krun_eff=krun_eff+1    !In case some replicates go extinct

c
c     Output:
c         

         OPEN(35,FILE='data_Ehux-EhV_PAPER_'//ck,Status='Unknown')
         write(35,17)'#Time(days)  Ehux(cell/mL)  EhuxI(cell/mL)  EhV(ce 
     &ll/mL) ALL_Ehux(cell/mL)  deaths(cell/mL)  MOI(-)  growth_rate(1/d
     &)  net_growth_rate(1/d)'

         DO i=1,itmax

              IF (points(i).ne.0d0) THEN
                 x1=time(i)/points(i)
                 x2=1d-3*EhuxF_t(i)/points(i)
                 x3=1d-3*EhuxI_t(i)/points(i)
                 x4=1d-3*EhV_t(i)/points(i)
                 x5=1d-3*Ehux_t(i)/points(i)
                 x6=1d-3*Ehux_deaths_t(i)/points(i)
                 x7=MOI_t(i)/points(i)
                 x8=mu_t(i)/points(i)
                 x9=net_mut(i)/points(i)

                 write(35,16) x1,x2,x3,x4,x5,x6,x7,x8,x9
                 
              END IF

           END DO
           CLOSE(35)
           OPEN(35,FILE='data_Ehux-EhV_PAPER_daily_'//ck,
     &          Status='Unknown')
         write(35,17)'#Time(days)  Ehux(cell/mL)  EhuxI(cell/mL)  EhV(ce 
     &ll/mL) ALL_Ehux(cell/mL)  deaths(cell/mL)  MOI(-)  growth_rate(1/d
     &)  net_growth_rate(1/d)'           
           t=1d0
           i=1
           it_fake=1d0
           DO WHILE (t.le.tmax)

              IF (mod(it_fake,steps+0d0).eq.0) i=i+1
              IF (points(i).ne.0d0) THEN
                 x1=time(i)/points(i)
                 x2=1d-3*EhuxF_t(i)/points(i)
                 x3=1d-3*EhuxI_t(i)/points(i)
                 x4=1d-3*EhV_t(i)/points(i)
                 x5=1d-3*Ehux_t(i)/points(i)
                 x6=1d-3*Ehux_deaths_t(i)/points(i)
                 x7=MOI_t(i)/points(i)
                 x8=mu_t(i)/points(i)
                 x9=net_mut(i)/points(i)

                 write(35,16) x1,x2,x3,x4,x5,x6,x7,x8,x9
                 
              END IF              
              t=t+dt
              it_fake=anint(it_fake+1d0)
           END DO
           CLOSE(35)           
      END DO

 16   FORMAT(9(E20.10,1x))
 17   FORMAT(A)
      end
