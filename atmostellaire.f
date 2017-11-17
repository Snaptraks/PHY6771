      program atmostellaire

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      real*4 time,dummy(2)
      dimension Tgris(100),Pgris(100) !
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/ !constante boltz,planck,vitesselum,masseelec
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/ !pi,mass H,c
      dimension deltaTsT(50,100),deltaHsH(50,100),condeq(50,100) !
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq  ! nu,w int,lambda
      common/couches/ tau(100),T(100),P(100),xKross(100),rho(100) !structure
      common/coeffscouchefreq/ xkappa(500,100),chi(500,100),sigma(100) !
      common/fluxall/ xJ(500,100),xH(500,100),xK(500,100),xHd(500,100) !
      common/planckdT/ dtau(500,100),plnk(500,100) !delta tau,planck (a chaque couche)
      common/correction/ deltaT(100),deltaH(100),deltaB(100)
      common/fluxmoyencouche/xJC(100),xBC(100),xHC(100) !
c
      open(21,file='atmo_6771_cooldau',status='old')
c
      nfreq=0.
      do i=1,500
         read(21,*,end=88) freq(i),poidsint(i),xlam(i)
         nfreq=nfreq+1.
      enddo
 88   continue
      close(21)
      
c 
c     Parametres modele
      Teff=10000.
      xlogg=4.0
      tau1=1.e-8
      tauND=1.e2
      ND=50
      sb=2.*(pi**5.)*(ek**4.)/(15.*(h**3.)*(c0**2.))
      Htot=sb*(Teff**4.)/(4.*pi)
c
c     Calcul de la structure grise
      call modelegris(tau1,tauND,Teff,xlogg,ND)
c     Garder en memoire la structure grise
      do id=1,ND
         Tgris(id)=T(id)
         Pgris(id)=P(id)
      enddo
C       write(*,*) 'structure grise done'
c
c     Calcul ETR pour modele gris
      call spectre(ND,xlogg)
C       write(*,*) 'ETR gris done'
c
c     Iteration pour correction de structure temp
      kble=0
      do k=1,50
c     Calcul condition equilibre
         do id=1,ND
            condeq(k,id)=0.
            do ij=1,nfreq
               condeq(k,id)=condeq(k,id)+poidsint(ij)*xkappa(ij,id)
     .              *(xJ(ij,id)-plnk(ij,id))
            enddo
            condeq(k,id)=condeq(k,id)/Htot
         enddo
c
C          write(*,*) '*************************************************'
C          write(*,*) 'iteration',k
         call correcT(xlogg,Htot,ND,sb)
         do id=1,ND
            deltaTsT(k,id)=deltaT(id)/(3.*T(id))
            deltaHsH(k,id)=deltaH(id)/Htot
         enddo
c         do id=1,ND
c            write(*,*) id,deltaT(id)/(3.*T(id)),deltaH(id)/Htot
c         enddo
c        Verifier convergence Temp
         do id=1,ND
            if (abs(deltaT(id)/(3.*T(id))).gt.1.e-3) goto 501
         enddo
         ktemp=k-kble
         kble=kble+1
C          write(*,*) 'Temperature convergee iteration',ktemp
c        Verifier convergence Flux
         do id=1,ND
            if (abs(deltaH(id)/Htot).gt.1.e-3) goto 501
         enddo
         kflux=k
C          write(*,*) 'Flux converge iteration',kflux
         goto 502
c        Recalculer structure
 501     call recalculstruct(ND,xlogg)
c        Calcul ETR pour modele corrige
         call spectre(ND,xlogg)
         if (k.eq.50) stop 'non convergence modele'
      enddo
c
 502  continue
c
c
      call etime(time,dummy)
      write(*,*) 'temps de calcul:',time
c
      end   !end program modeleHpure
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine eqetat(P,T)
      implicit real*8 (a-h,o-z)
      dimension gH(16),chiN(16)              ! gn et chin pour HI
      data chiH1,chiHm/157896.0,8763.0/ !energie ionisation
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/
      data pi,xmH/3.141592654,1.67262310d-24/
      data uH2,uHm,tol/1.,1.,1.d-7/
      common/pops/ xN(4),xNH1(16),rhod   ! xN =(Ne,NH1,NH2,NH-),xNH1=population niv energie,densite

c
c     Calcul xNtot
      xNtot = P/(ek*T)
c
c     Calcul de gn et chiN et uH1 n=1,5
      uH1 = 0.
      do n = 1,16
         gH(n) = 2.*(n**2.)
         chiN(n) = chiH1*(1.-1./(n**2.))  !la constante de boltzman est incluse dans le chi
         uH1 = uH1 + gH(n)*exp(-chiN(n)/T) 
      enddo
c
c     Calcul de phi1 et phi2
      A = (2.*pi*xme*ek*T/(h**2.))**(3./2.)
      phi1 = ((2.*A/uH1)*exp(-chiH1/T))**(-1.)
      phi2 = ((2.*A*uH1)*exp(-chiHm/T))**(-1.)
c
c     Calcul de Ne par Newton-Rawphson (boucle iterative)
      xN(1) = ((1.+phi1*xNtot)**(1./2.)-1.)/phi1 !on pose une valeur initiale de Ne0
c
      do i = 1,100
         TPH = 1.-(xN(1)**2.)*phi1*phi2
         BTH = 1.+xN(1)*phi1+(xN(1)**2.)*phi1*phi2
         F = xN(1) - (xNtot-xN(1))*(TPH/BTH)
         dFdNe = 1.+(TPH/BTH)+(xNtot-xN(1))*(2.*xN(1)*phi1*phi2*BTH  ! calcul dF/dNe
     .           +(phi1+2.*xN(1)*phi1*phi2)*TPH)/(BTH**2.)
         xN(1) = xN(1)-F/dFdNe       ! on corrige Ne
         if (abs(F/(dFdNe*xN(1))).lt.tol) goto 201 ! si F<tolerance,sort boucle
         if (i.eq.100) stop 'non convergence eqetat'
      enddo
c
c     Calcul des autres populations
c
 201  xN(3) = xN(1)/TPH
      xN(4) = xN(3)*phi1*(xN(1)**2.)*phi2
      xN(2) = xNtot - xN(1) - xN(3) - xN(4)
c
c     Calcul des populations des niv energie HI
c
      do n = 1,16
         xNH1(n) = (xN(2)*gH(n)/uH1)*exp(-chiN(n)/T)
      enddo
c
c     Calcul de la densite de masse
c      
      rhod = (xNtot-xN(1))*xmH   !nombre de noyaux H* masse H neutre
c      
      return
c
      end   !end subroutine eqetat
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine opac(P,T)
      implicit real*8 (a-h,o-z)
      real rT       ! pression gazeuse,temperature,racineT
      real nu0(16)  !freq coupure 
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/
      data pi,xmH/3.141592654,1.67262310d-24/
      data A,B,C,alphadiff/2.815d29,3.29d15,3.69d8,6.6516d-25/  !constantes pour les alpha
      common/pops/ xN(4),xNH1(16),rhod   ! xN =(Ne,NH1,NH2,NH-),xNH1=population niv energie,densite
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq 
      common/coeff/ xkappa(500),chi(500),sigma  !/rho
c
c     Calcul des populations
      call eqetat(P,T)
c
c     Calcul de alphabf,alphaff,kappa,sigma,chi
      rT = T**(1./2.)
      sigma=xN(1)*alphadiff/rhod
c      
      do j=1,nfreq
         xkappa(j)=xN(1)*xN(3)*C/((freq(j)**3.)*rT)
         do i=1,16
            nu0(i)=B/(i**2.)
            if (freq(j).ge.nu0(i)) then
               xkappa(j)=xkappa(j)+xNH1(i)*A/((i**5.)*(freq(j)**3.))
            endif
         enddo
         xkappa(j)=xkappa(j)*(1.-exp(-1.*h*freq(j)/(ek*T)))/rhod
         chi(j)=xkappa(j)+sigma
      enddo
c
      return
c         
      end !opac
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ross(P,T,opacR)
      implicit real*8 (a-h,o-z)
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/   ! 
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq ! nu,w int,lambda
      common/coeff/ xkappa(500),chi(500),sigma  !sections efficaces/rho
c
c     Calcul de l'opacite
      call opac(P,T)
c
c     Calcul des integrales et opacR
      C=h/(ek*T)
      top=0.
      bot=0.
      do j=1,nfreq
         B=(2./ek)*((h/(c0*T))**2.)*(freq(j)**4.)*exp(freq(j)*C)
     .     *(exp(freq(j)*C)-1)**(-2.)
         top=top+poidsint(j)*(1./chi(j))*B
         bot=bot+poidsint(j)*B
      enddo
      opacR=bot/top
c      
      return
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine modelegris(tau1,tauND,Teff,xlogg,ND)
      implicit real*8 (a-h,o-z)
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/ !
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/    !
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq !
      common/couches/ tau(100),T(100),P(100),xKross(100),rho(100)
c
      tau(1)=tau1
      tau(ND)=tauND
      deltaTau=(log10(tauND)-log10(tau1))/(ND-1)
c
c     Calcul echelle tau
      do i=2,ND-1
         tau(i)=10.**(log10(tau(i-1))+deltaTau)
      enddo
c
c     Calcul structure Temp
      do i=1,ND
         q=0.7104-0.1331*exp(-3.4488*tau(i))
         T(i)=((3./4.)*(Teff**4.)*(tau(i)+q))**(1./4.)
      enddo
c
c     Calcul structure Press
      tol=10.**(-6.)
      P(1)=1.
      do j=1,30
         call ross(P(1),T(1),xKross(1))
         xNewP=(10.**xlogg)*tau(1)/xKross(1)
         deltaP=abs((xNewP-P(1))/xNewP)
         P(1)=xNewP
         if (deltaP.lt.tol) goto 202
         if (j.eq.30) stop 'non convergence P(1)'
      enddo
 202  continue
C       write(*,*) 'P(1) done'
c
      tol=10**(-6.)
      do i=2,ND
         P(i)=P(i-1)+(10.**xlogg)*(tau(i)-tau(i-1))/xKross(i-1)
         do j=1,100
            call ross(P(i),T(i),xKross(i))
            xNewP=P(i-1)+2.*(10.**xlogg)*(tau(i)-tau(i-1))
     .             /(xKross(i-1)+xKross(i))
            deltaP=abs((xNewP-P(i))/P(i))
            P(i)=xNewP
            if (deltaP.lt.tol) goto 203
            if (j.eq.100) stop 'non convergence P(i)'
         enddo
 203     continue
      enddo
C       write(*,*) 'structure P done'
c
      return
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine matinv(a,n,nr)
      implicit real*8 (a-h,o-z)
      dimension a(nr,nr)
      data one/1.0/
c     ******************************************************************                                    
c     Matrix inversion of a by LU reduction.                                                                
c     LU decomposition , starting with L.                                                                   
c     ****************************************************************** 
      do 5 i = 2, n
         im1 = i - 1
         do 2 j = 1, im1
            jm1 = j - 1
            div=a(j,j)
            sum = 0.0e0
            if(jm1 .lt. 1 ) go to 2
            do 1 l = 1, jm1
    1       sum = sum + a(i,l)*a(l,j)
    2       a(i,j)=(a(i,j)-sum)/div
         do 4 j = i,n
            sum = 0.0e0
            do 3 l = 1, im1
    3       sum = sum + a(i,l)*a(l,j)
            a(i,j) = a(i,j) - sum
    4    continue
    5 continue
c     ******************************************************************                                    
c     Inversion of L                                                                                        
c     ******************************************************************  
      do 13 ii = 2, n
         i = n + 2 - ii
         im1 = i - 1
         if(im1.lt.1) go to 13
         do 12 jj = 1, im1
            j = i - jj
            jp1 = j + 1
            sum = 0.0e0
            if(jp1.gt.im1) go to 12
            do 11 k = jp1, im1
   11       sum = sum + a(i,k)*a(k,j)
   12       a(i,j) = - a(i,j) - sum
   13 continue
c     ******************************************************************                                    
c     U inversion                                                                                           
c     ******************************************************************                                    
      do 17 ii = 1, n
         i = n + 1 - ii
         div=a(i,i)
         ip1 = i + 1
         if(ip1.gt.n) go to 17
         do 16 jj = ip1, n
            j = n + ip1 - jj
            sum = 0.0e0
            do 15 k = ip1, j
   15       sum = sum + a(i,k)*a(k,j)
            a(i,j)=-sum/div
   16    continue
   17 a(i,i) = one/a(i,i)
c     ******************************************************************                                    
c     Multiplication of U inverse and L inverse                                                             
c     ******************************************************************                                    
      do 24 i = 1, n
         do 23 j = 1, n
           k0= max0(i,j)
            if(k0.eq.j) go to 22
            sum = 0.0d0
   20       do 21 k = k0, n
   21       sum = sum + a(i,k)*a(k,j)
            go to 23
   22       sum = a(i,k0)
            if(k0.eq.n) go to 23
            k0 = k0 + 1
            go to 20
   23    a(i,j) = sum
   24 continue
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine feautrier(ND)
      implicit real*8 (a-h,o-z)
      parameter (M=3)
      dimension xmu(3),w(3)               !angles
      dimension A(3),B(3,3),C(3),xL(3),T(3,3),X(3)
      dimension u(3,100),R(3,100),D(3,3,100) !matrices
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/ !
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/ !
      data xmu/0.8872983346,0.5,0.1127016654/ !
      data w/0.27777777,0.4444444,0.27777777/ !
      common/flux/ xJ(100),xH(100),xK(100)
      common/plnkdTlamnu/ dtau(100),plnk(100),xlamnu(100) !
c
c
c     Couche 1
      do j=1,M
         C(j)=-xmu(j)/dtau(2)
         xL(j)=(-0.5)*dtau(2)*xlamnu(1)*plnk(1)/xmu(j)
         do k=1,M
            B(j,k)=0.5*dtau(2)*(1.-xlamnu(1))*w(k)/xmu(j)
         enddo
         B(j,j)=B(j,j)+C(j)-1.-0.5*dtau(2)/xmu(j)
      enddo
c     inverser B
      call matinv(B,M,M)
c     calcul de D(1)
      do j=1,M
         R(j,1)=0.
         do k=1,M
            D(j,k,1)=B(j,k)*C(k)
            R(j,1)=R(j,1)+B(j,k)*xL(k)
         enddo
      enddo
c     Couches 2,ND-1
      do id=2,ND-1
         do j=1,M
            A(j)=(-2.)*(xmu(j)**2.)/(dtau(id)*(dtau(id)+dtau(id+1)))
            C(j)=(-2.)*(xmu(j)**2.)/(dtau(id+1)*(dtau(id)+dtau(id+1)))
            xL(j)=(-1.)*xlamnu(id)*plnk(id)
            do k=1,M
               B(j,k)=(1.-xlamnu(id))*w(k)
            enddo
            B(j,j)=B(j,j)+A(j)+C(j)-1.
         enddo
         do j=1,M
            do k=1,M
               T(j,k)=B(j,k)-A(j)*D(j,k,id-1)
            enddo
            X(j)=A(j)*R(j,id-1)+xL(j)
         enddo
         call matinv(T,M,M)
         do j=1,M
            R(j,id)=0.
            do k=1,M
               D(j,k,id)=T(j,k)*C(k)
               R(j,id)=R(j,id)+T(j,k)*X(k)
            enddo
         enddo
      enddo
c     Couche ND
      do j=1,M
         u(j,ND)=plnk(ND)
      enddo
c     Subtitution inverse
      do id=ND-1,1,-1
         do j=1,M
            u(j,id)=0.
            do k=1,M
               u(j,id)=u(j,id)+D(j,k,id)*u(k,id+1)
            enddo
            u(j,id)=u(j,id)+R(j,id)
         enddo
      enddo
c     Calcul J,H,K
      do id=1,ND
         xJ(id)=0.
         xK(id)=0.
         do j=1,M
            xJ(id)=xJ(id)+w(j)*u(j,id)
            xK(id)=xk(id)+w(j)*(xmu(j)**2.)*u(j,id)
         enddo
      enddo
      do id=2,ND
         xH(id)=(xK(id)-xK(id-1))/dtau(id)
      enddo
      xH(1)=0.
      do j=1,M
         xH(1)=xH(1)+w(j)*xmu(j)*u(j,1)
      enddo
c
      return
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine spectre(ND,xlogg)
      implicit real*8 (a-h,o-z)
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/ !
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/ !
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq  ! nu,w int,lambda
      common/flux/ xJnu(100),xHnu(100),xKnu(100) !
      common/fluxall/ xJ(500,100),xH(500,100),xK(500,100),xHd(500,100) !
      common/coeff/ xkappad(500),chid(500),sigmad  !sections efficaces
      common/coeffscouchefreq/ xkappa(500,100),chi(500,100),sigma(100) !
      common/planckdT/ dtau(500,100),plnk(500,100)         !delta tau,planck (a chaque couche)
      common/plnkdTlamnu/ dtaunu(100),plnknu(100),xlamnu(100) !
      common/couches/ tau(100),T(100),P(100),xKross(100),rho(100) !
      common/pops/ xN(4),xNH1(16),rhod   ! xN =(Ne,NH1,NH2,NH-),xNH1=population niv energie,densite
c
c     Calcul des opac chaque couche
      do id=1,ND
         call opac(P(id),T(id))
         do ij=1,nfreq
            xkappa(ij,id)=xkappad(ij)
            chi(ij,id)=chid(ij)
            sigma(id)=sigmad
            rho(id)=rhod
         enddo
      enddo
c
c     Boucle sur les frequences
      do ij=1,nfreq
c        Calcul delta tau
         dtaunu(1)=chi(ij,1)*P(1)/(10.**xlogg)
         do id=2,ND
            dtaunu(id)=0.5*(chi(ij,id)+chi(ij,id-1))
     .           *(P(id)-P(id-1))/(10.**xlogg)
         enddo
         do id=1,ND
c           calcul de Bnu
            plnknu(id)=2.*h*(freq(ij)**3.)/
     .           ((c0**2.)*(exp(h*freq(ij)/(ek*T(id)))-1.))
c           calcul de lambdanu
            xlamnu(id)=xkappa(ij,id)/chi(ij,id)
         enddo
c        Resolution ETR
         call feautrier(ND)
c         write(*,*) 'TR done point freq',ij
c        Sauvegarde dTau,j,h,k pour tout couche et freq
         do id=1,ND
            dtau(ij,id)=dtaunu(id)
            plnk(ij,id)=plnknu(id)
            xJ(ij,id)=xJnu(id)
            xK(ij,id)=xKnu(id)
            xH(ij,id)=xHnu(id)
         enddo
c        Calcul de H aux couches a partir des demi-couches
         xHd(ij,1)=xH(ij,1)
         xHd(ij,ND)=xH(ij,ND)
         do id=2,ND-1
            xHd(ij,id)=(xH(ij,id)+xH(ij,id+1))/2.
         enddo
      enddo
c
      return
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine correcT(xlogg,Htot,ND,sb)
      implicit real*8 (a-h,o-z)
      dimension xkapJ(100),xkapP(100),chiF(100)    !
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/ !
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/ !
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq  ! nu,w int,lambda
      common/fluxall/ xJ(500,100),xH(500,100),xK(500,100),xHd(500,100) !
      common/coeffscouchefreq/ xkappa(500,100),chi(500,100),sigma(100) !
      common/planckdT/ dtau(500,100),plnk(500,100)         !delta tau,planck (a chaque couche)
      common/correction/ deltaT(100),deltaH(100),deltaB(100)
      common/couches/ tau(100),T(100),P(100),xKross(100),rho(100) !
      common/fluxmoyencouche/xJC(100),xBC(100),xHC(100) !
c

c     Calcul deltaH
      do id=1,ND
         deltaH(id)=Htot
         do ij=1,nfreq
            deltaH(id)=deltaH(id)-(poidsint(ij)*xHd(ij,id))
         enddo
      enddo
c
      do id=1,ND
c        Calcul de J,B,H couches
         xJC(id)=0.
         xBC(id)=0.
         xHC(id)=0.
         do ij=1,nfreq
            xJC(id)=xJC(id)+poidsint(ij)*xJ(ij,id)
            xBC(id)=xBC(id)+poidsint(ij)*plnk(ij,id)
            xHC(id)=xHC(id)+poidsint(ij)*xHd(ij,id)
         enddo
c
c        Calcul kappaJ,kappaP,chiF
         xkapJ(id)=0.
         xkapP(id)=0.
         chiF(id)=0.
         do ij=1,nfreq
            xkapJ(id)=xkapJ(id)+poidsint(ij)*xkappa(ij,id)*xJ(ij,id)
            xkapP(id)=xkapP(id)+poidsint(ij)*xkappa(ij,id)*plnk(ij,id)
            chiF(id)=chiF(id)+poidsint(ij)*chi(ij,id)*xHd(ij,id)
         enddo
         xkapJ(id)=xkapJ(id)*rho(id)/xJC(id)
         xkapP(id)=xkapP(id)*rho(id)/xBC(id)
         chiF(id)=chiF(id)*rho(id)/xHC(id)
      enddo
c
c     Calcul deltaB
      do id=1,ND
         deltaB(id)=0.
         if (id.eq.1) goto 300
         do iid=1,id-1
            deltaB(id)=deltaB(id)+(3./(2.*(10.**xlogg)))
     .           *((chiF(iid)*deltaH(iid)/rho(iid))
     .           +(chiF(iid+1)*deltaH(iid+1)/rho(iid+1)))
     .           *(P(iid+1)-P(iid))
         enddo
 300     deltaB(id)=(xkapJ(id)/xkapP(id))
     .        *(xJC(id)+deltaB(id)+2.*deltaH(1))-xBC(id)
c
c        Calcul deltaT
         deltaT(id)=(pi/(4.*sb*(T(id)**3.)))*deltaB(id)
      enddo
c
      return
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine recalculstruct(ND,xlogg)
      implicit real*8 (a-h,o-z)
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/  !
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/ !   
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq !
      common/couches/ tau(100),T(100),P(100),xKross(100),rho(100) !
      common/correction/ deltaT(100),deltaH(100),deltaB(100)
c
c     Changement temp      
      do id=1,ND
         T(id)=T(id)+(deltaT(id)/3.)
      enddo      
c
c     Calcul structure Press
      tol=1.e-6
      P(1)=1.
      do j=1,30
         call ross(P(1),T(1),xKross(1))
         xNewP=(10.**xlogg)*tau(1)/xKross(1)
         deltaP=abs((xNewP-P(1))/xNewP)
         P(1)=xNewP
         if (deltaP.lt.tol) goto 602
         if (j.eq.30) stop 'non convergence P(1)'
      enddo
 602  continue
c
      do id=2,ND
         P(id)=P(id-1)+(10.**xlogg)*(tau(id)-tau(id-1))/xKross(id-1)
         do j=1,30
            call ross(P(id),T(id),xKross(id))
            xNewP=P(id-1)+2.*(10.**xlogg)*(tau(id)-tau(id-1))
     .             /(xKross(id-1)+xKross(id))
            deltaP=abs((xNewP-P(id))/P(id))
            P(id)=xNewP
            if (deltaP.lt.tol) goto 603
            if (j.eq.30) stop 'non convergence P(i)'
         enddo
 603     continue
      enddo
c
      return
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function fctpart(Z, N, T)
C       Calcule la fonction de partition pour un atome de numéro atomique Z
C       de nombre d'électron N et à une température T.
      implicit real*8 (a-h,o-z)
      integer Z, N
      
      fctpart = 1.
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function ionlevel(Z, n)
C       Retourne les énerdies d'ionisations de l'atome de numéro atomique Z
C       pour le niveau d'ionisation n.
      integer Z, n
      
      ionlevel = 10. ! eV
      
      return
      end