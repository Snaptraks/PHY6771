      program atmostellaire

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      real*4 time,dummy(2)
      dimension Tgris(100),Pgris(100) !
      dimension Az(6)
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/ !constante boltz,planck,vitesselum,masseelec
      data pi,xmH,c0/3.141592654,1.67262310d-24,2.99792458d10/ !pi,mass H,c
      dimension deltaTsT(50,100),deltaHsH(50,100),condeq(50,100) !
      common/abond/ Az ! abondances de O,Ne,Na,Mg,Al,Si
      common/ptsfreq/ freq(500),poidsint(500),xlam(500),nfreq  ! nu,w int,lambda
      common/couches/ tau(100),T(100),P(100),xKross(100),rho(100) !structure
      common/coeffscouchefreq/ xkappa(500,100),chi(500,100),sigma(100) !
      common/fluxall/ xJ(500,100),xH(500,100),xK(500,100),xHd(500,100) !
      common/planckdT/ dtau(500,100),plnk(500,100) !delta tau,planck (a chaque couche)
      common/correction/ deltaT(100),deltaH(100),deltaB(100)
      common/fluxmoyencouche/xJC(100),xBC(100),xHC(100) !
      
      common/pops/ xNe,xNtot,xNz(6,3),rhod   ! xNive(6,3,?) tableau niv energie.. ?? combien de niveaux
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
      Az(1)=10.**-0.00 ! O
      Az(2)=10.**-0.00 ! Ne
      Az(3)=10.**-2.00 ! Na
      Az(4)=10.**-2.00 ! Mg
      Az(5)=10.**-2.50 ! Al
      Az(6)=10.**-3.00 ! Si
      
      sb=2.*(pi**5.)*(ek**4.)/(15.*(h**3.)*(c0**2.))
      Htot=sb*(Teff**4.)/(4.*pi)
      
      ! Lire les fichiers topbase
      call initiateEnergyLevels()
      
      open(10, file='out_test_6.txt')
      
      do i = 70,1200
         call eqetat(1d4, i*50d0)
         write(10, '(4e13.4)') i*50d0, xNz(6,:)/xNtot
      enddo
      close(10)
      
      stop 'I CANCELLED IT'
c
c     Calcul de la structure grise
      call modelegris(tau1,tauND,Teff,xlogg,ND)
c     Garder en memoire la structure grise
      do id=1,ND
         Tgris(id)=T(id)
         Pgris(id)=P(id)
      enddo
      write(*,*) 'structure grise done'
      stop
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
      end   !end program atmostellaire
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine eqetat(P,T)
      implicit real*8 (a-h,o-z)
      integer Z, xNelec
      real*8 fctpart, ionlevel
      dimension AZ(6),xNelec(6,3) !
      dimension phi1(6),phi2(6) !
      dimension chiz(6,2),U(6,3)  ! (divisé par k en electron volt)
      dimension TP(6),BT(6) ! top et bottom
      dimension xmass(6),Z(6)
      data ek,h,xme/1.38065d-16,6.6260755d-27,9.1093897d-28/
      data pi/3.141592654/
      ! masses des atomes, selon Google "mass <element> atom in grams"
      data xmass/2.6566962d-23, ! O
     .           3.3509177d-23, ! Ne
     .           3.8175407d-23, ! Na
     .           4.0359398d-23, ! Mg
     .           4.4803895d-23, ! Al
     .           4.6637066d-23/ ! Si
      data tol/1.d-7/
      data Z/8,10,11,12,13,14/ !O,Ne,Na,Mg,Al,Si
      common/abond/ Az ! abondances de O,Ne,Na,Mg,Al,Si
      common/pops/ xNe,xNtot,xNz(6,3),rhod   ! xNive(6,3,?) tableau niv energie.. ?? combien de niveaux
      
c
c     Ecrire xNelec, nombre d'electron pour chaque espece
      do i = 1,6
         do k = 1,3 
            xNelec(i,k)=Z(i)-(k-1)! Atome neutre a tous ses electrons
         enddo
      enddo
c
c     Appeler subroutine fct partition
      do i = 1,6
         do k=1,3
            U(i,k) = fctpart(Z(i),xNelec(i,k),T)
            chiz(i,k) = ionlevel(Z(i), xNelec(i, k))
C             print*, Z(i), xNelec(i, k), U(i, k), chiz(i, k)
         enddo
      enddo
c     
c     Calcul xNtot
      xNtot = P/(ek*T)
c
c     Calcul de phi1 et phi2
      A = (2.*pi*xme*ek*T/(h**2.))**(3./2.)
      do i=1,6
         phi1(i) = ((2.*A*U(i,2)/U(i,1))*exp(-chiz(i,1)/T))**(-1.)
         phi2(i) = (2.*A*U(i,3)/U(i,2))*exp(-chiz(i,2)/T)
      enddo
c
c     Calcul de Ne par Newton-Rawphson (boucle iterative)
      !on pose une valeur initiale de Ne0 comme si O pur
      xNe = ((1.+phi1(1)*xNtot)**(1./2.)-1.)/phi1(1)
c
      do j = 1,100
         sum1 = 0.
         sum2 = 0.
         do i = 1,6
            TP(i) = xNe+2.*phi2(i)
            BT(i) = phi2(i)+xNe+(xNe**2.)*phi1(i)
            sum1 = sum1+Az(i)*TP(i)/BT(i)
            sum2 = sum2+Az(i)*(BT(i)-TP(i)*(1+2.*xNe*phi1(i)))
     .           /(BT(i)**2.)
         enddo
         F = xNe-(xNtot-xNe)*sum1
         dFdNe = 1.+sum1-(xNtot-xNe)*sum2
         
         xNe = xNe-F/dFdNe       ! on corrige Ne
         
         if (abs(F/(dFdNe*xNe)).lt.tol) goto 201 ! si F<tolerance,sort boucle
         if (j.eq.10) stop 'non convergence eqetat'
      enddo
c
c     Calcul des autres populations
c
 201  continue
c     
      do i = 1,6
         xNz(i,2) = (xNtot-xNe)*Az(i)*xNe/
     .              (phi2(i)+xNe+(xNe**2.)*phi1(i))
         xNz(i,1) = xNe*xNz(i,2)*phi1(i)
         xNz(i,3) = xNz(i,2)*phi2(i)/xNe
      enddo
c
c     Calcul des populations des niv energie HI
c
c     ??? combien de niveaux on met par espece
c
c     Calcul de la densite de masse
c      
      rhod = 0.
      do i = 1,6
         do k = 1,3
            rhod = rhod+xNz(i,k)*xmass(i)
         enddo
      enddo
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
C       print*, bot, top
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
      write(*,*) 'P(1) done'
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
      function ionlevel(Z, n) result(chi)
C       Retourne les energies d'ionisations de l'atome de numero atomique Z
C       avec n electrons (ionise Z-n fois).
      integer, intent(in) :: Z, n
      real*8 chi, keV
      integer i
      real*8, dimension(6,3) :: energyLevel
      integer, dimension(6)  :: ielem
      
      ! CRC Handbook of Chemestry and Physics (eV)
      data energyLevel(1,:)/13.61805, 35.1211,  54.9355 / ! O
      data energyLevel(2,:)/21.56454, 40.96296, 63.45   / ! Ne
      data energyLevel(3,:)/5.139076, 47.2864,  71.6200 / ! Na
      data energyLevel(4,:)/7.646235, 15.03527, 80.1437 / ! Mg
      data energyLevel(5,:)/5.985768, 18.82855, 28.44765/ ! Al
      data energyLevel(6,:)/8.15168,  16.34584, 33.49302/ ! Si
      
      data ielem/8, 10, 11, 12, 13, 14/
      data keV  / 8.6173303D-5 /!eV/k  | Constante de Boltzman
      
      do i=1,6
         if (ielem(i) == Z) exit
      enddo
      
      chi = energyLevel(i, Z-n+1)/keV
      
      end function ionlevel
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ! Lit les fichiers de topbase et entre les valeurs dans un
      ! block common qui peut être réutilisé plus tard.
      ! nbLevels = nombre de niveaux d'énerie pour chaque ions
      ! energy   = énergie des différents niveaux pour chaque ions
      ! g        = poinds statistique pour calculer les fonctions de partitions
      SUBROUTINE initiateEnergyLevels()
      
      implicit none
      integer,dimension(20,2)  :: ionList
      integer,dimension(20)    :: nbLevels
      real*8,dimension(20,500) :: energy, g
      integer        :: nbIons, atomicNumber, electronNumber
      integer        :: i, j, k, NZread, NEread, nbl, io
      character*6    :: filename, tmpChar
      
      common/ions/ ionList, nbLevels, energy, g, nbIons
      
      !List of ions
      ionList(1,:) =  [8 ,6 ] !OIII
      ionList(2,:) =  [8 ,7 ] !OII
      ionList(3,:) =  [8 ,8 ] !OI
      ionList(4,:) =  [10,8 ] !NeIII
      ionList(5,:) =  [10,9 ] !NeII
      ionList(6,:) =  [10,10] !NeI
      ionList(7,:) =  [11,9 ] !NaIII
      ionList(8,:) =  [11,10] !NaII
      ionList(9,:) =  [11,11] !NaI
      ionList(10,:) = [12,10] !MgIII
      ionList(11,:) = [12,11] !MgII
      ionList(12,:) = [12,12] !MgI
      ionList(13,:) = [13,11] !AlIII
      ionList(14,:) = [13,12] !AlII
      ionList(15,:) = [13,13] !AlI
      ionList(16,:) = [14,12] !SiIII
      ionList(17,:) = [14,13] !SiII
      ionList(18,:) = [14,14] !SiI
      nbIons = 18
      
      !Ouvrir et lire les fichiers topbase
      do i=1,nbIons
         atomicNumber = ionList(i,1)
         electronNumber = ionList(i,2)
         
         if (atomicNumber < 10) then !Trouver le nom du fichier, ZXXEXX
            write (filename,'(A2,I1)') 'Z0',atomicNumber
         else
            write (filename,'(A1,I2)') 'Z',atomicNumber
         end if
         if (electronNumber < 10) then
            write (tmpChar,'(A2,I1)') 'E0',electronNumber
         else
            write (tmpChar,'(A1,I2)') 'E',electronNumber
         end if
         
         filename = trim(filename)//trim(tmpChar) !Ouvrir le fichier
         open(211,file='topbase/'//filename,status='old')
         read(211,*) 
         read(211,*) 
         read(211,*) 
         
         do j=1,500 !Lire le fichier
            read(211,'(8X,2I3,37X,2E13.5)',IOSTAT=io) NZread, NEread,
     $                                            energy(i,j), g(i,j)
            if (io < 0) EXIT !End of file
            !Erreurs possibles
            if (NZread /= atomicNumber .or. NEread /=  electronNumber) 
     $             STOP 'Error during topbase files reading.
     $                   Atomic/electron number mismatch.'
     
            if (j>1 .and. energy(i,j) < energy(i,j-1)) 
     $             STOP 'Error during topbase files reading.
     $                   Energy levels in wrong order.'

         end do
         nbLevels(i) = j-1
      
      end do
      
      END SUBROUTINE initiateEnergyLevels
      
      
      ! Calcule la fonction de partition pour une température T (real*8),
      ! un numéro atomique atomicNumber (integer) et un nombre d'électron
      ! electronNumber (integer)
      ! La sous-routine initiateEnergyLevels doit avoir été appelée une
      ! fois avant.
      FUNCTION fctpart(atomicNumber, electronNumber, T) result(U)
      
      implicit none
      integer, intent(in)      :: atomicNumber, electronNumber
      real*8, intent(in)       :: T
      integer,dimension(20,2)  :: ionList
      integer,dimension(20)    :: nbLevels
      real*8,dimension(20,500) :: energy, g
      real*8                   :: Ryd, keV, U, partitionHydrogene
      integer                  :: i, j, nbIons
      data Ryd  / 13.605693D0  /!eV    | Rydberg
      data keV  / 8.6173303D-5 /!eV/k  | Constante de Boltzman
      
      common/ions/ ionList, nbLevels, energy, g, nbIons
      
      !Hydrogène
      U = 0.d0
      if (atomicNumber==1 .and. electronNumber==1) then 
         U = partitionHydrogene(T,16)
         
      !Les autres ions 
      else 
         U = 0.d0                                          
         do i=1,nbIons+1 !Trouver l'index dans les tableaux
            if (ionList(i,1) == atomicNumber .and. 
     $          ionList(i,2) == electronNumber) EXIT
         end do
         
         if (i==nbIons+1) STOP 'Ion does not exit'
         
         do j=1,nbLevels(i) !Calculer U
            U = U + g(i,j)* EXP(-energy(i,j)*Ryd/keV/T )
         end do
      end if

      END FUNCTION fctpart
      
      !Fonction de partition de l'hydrogène
      !Notes de Gilles Fontaine, éq. 1.19, page 16
      FUNCTION partitionHydrogene(T, nlimit) result(U)
         
         implicit none
         integer, intent(in)  :: nlimit
         real*8, intent(in)   :: T
         integer              :: i, j
         real*8               :: U, Ryd, keV
         data Ryd  / 13.605693D0  /!eV    | Rydberg
         data keV  / 8.6173303D-5 /!eV/k  | Constante de Boltzman

         U = 0.d0
         do i=1,nlimit
            U = U + 2.d0*(DBLE(i)**2) * EXP(-Ryd/keV/T * 
     $          (1.d0-1.d0/(i**2)))
         end do

      END FUNCTION partitionHydrogene
      
      
      
