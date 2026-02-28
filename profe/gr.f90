subroutine calcul_gr(nmolecules,natoms,r,costat)
      implicit double precision(a-h,o-z)
      include 'exercicishake.dim'
      dimension r(3,nmax,nmaxmol)
      
      parameter (nbins = 200)
      dimension hist(nbins)
      
      pi = 4.d0*datan(1.d0)
      
c     Inicialitzar l histograma
      do i = 1, nbins
         hist(i) = 0.d0
      end do
      
c     El rang maxim es la meitat de la caixa per evitar auto-interaccions
      rmax = costat/2.d0
      del_r = rmax/dfloat(nbins)
      
c     Bucle sobre tots els parells d atoms de molecules DIFERENTS
      do ic = 1,nmolecules-1
         do is = 1,natoms
            do jc = ic+1,nmolecules
               do js = 1,natoms
                  
c                 Calcul de distancia amb minima imatge
                  rr2 = 0.d0
                  do l = 1,3
                     rijl = r(l,js,jc) - r(l,is,ic)
                     rijl = rijl - costat*dnint(rijl/costat)
                     rr2 = rr2 + rijl*rijl
                  end do
                  dist = dsqrt(rr2)

c                 Clasificar dins l histograma si esta en el rang
                  if (dist.lt.rmax) then
                     ibin = dint(dist/del_r) + 1
                     if (ibin.le.nbins) then
                        hist(ibin) = hist(ibin) + 2.d0 ! Parell (i,j) i (j,i)
                     end if
                  end if

               end do
            end do
         end do
      end do

c     Normalitzacio y escriptura de resultats
      vol = costat**3
c     Densitat efectiva para parells intermoleculars
      rho_inter = dfloat((nmolecules-1)*natoms)/vol
      total_atoms = dfloat(nmolecules*natoms)

      open(60,file='gr.dat',status='unknown')
      do i = 1,nbins
         r_bin = (dfloat(i)-0.5d0)*del_r
c        Volum de la closca esferica exacta
         vol_shell = (4.d0/3.d0)*pi*((r_bin+0.5d0*del_r)**3 - 
     &               (r_bin-0.5d0*del_r)**3)
         
         ideal_count = total_atoms * rho_inter * vol_shell
         g_val = hist(i) / ideal_count
         
         write(60,*) r_bin, g_val
      end do
      close(60)

      return
      end
