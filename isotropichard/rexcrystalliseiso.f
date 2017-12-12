      subroutine fisotropic(sy,dsy,depl,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
c
      implicit real*8(a-h,o-z)
      logical checkrx
      dimension props(nprops),statev(nstatv),tempstatev(nstatv),
     & xi(2),xj(2),r(2),drdx(2,2),fxinfo(3),reps(2),dxdr(2,2),
     & fxnvec(5),xjupd(2)
     
      parameter(zero=0.d0,half=0.5d0,one=1.d0,two=2.d0,
     & toler=1.d-4,x10=one,x20=1.d-10,fxn0=1.d-4,ratelim=1.d-8)
c
c elastic properties:
      emu0 = props(1)
      ed0 = props(2) 
      et0 = props(3) 
      enu = props(4)
c reference stress values
      siga = props(5)
      sig0 = props(6)
c scaling function
      a0e = props(7)
      rate0 = props(8)
      qe = props(9)
      pe = props(10)
c
      rate = depl/dtime
      if(rate.lt.ratelim)then
      rate = ratelim
      endif
c
      if(temp.gt.et0)then
      emu = emu0 - ed0/(dexp(et0/temp)-one)
      sfe0 = temp/(a0e*emu)
      else
      emu = emu0
      sfe0 = one/a0e
      endif
      emusf = emu/emu0
c 
      sfel = dlog(rate0/rate)*sfe0
      sfe = dabs(one-sfel**(one/qe))**(one/pe)
c isv shift
      nrrx = (nstatv-3)/4 
      if(statev(7).gt.(0.999d0))then
      ixvf=0
      do while(ixvf.lt.nrrx)
      lstskip = 4*(ixvf)+3
      statev(lstskip+1)=statev(lstskip+5)
      statev(lstskip+2)=statev(lstskip+6)
      statev(lstskip+3)=statev(lstskip+7)
      statev(lstskip+4)=statev(lstskip+8)
      enddo
      statev(lstskip+5)=zero
      statev(lstskip+6)=zero
      statev(lstskip+7)=zero
      statev(lstskip+8)=zero
      endif
c    
      fxc = one
      fxcp = one
      fxcr = zero
      dfxcde = zero
      x1eq = zero
      x2eq = zero
      dx1eqde = zero
      plastic = zero
c
      ixvf = 1
      checkrx = .true.
      do while((ixvf.lt.nrrx).and.(checkrx))
      lstskip = 4*(ixvf-1)+3
      xeplp = statev(lstskip+1)
      x1p = max(statev(lstskip+2),x10)
      x2p = max(statev(lstskip+3),x20)
      fxnp = max(statev(lstskip+4),fxn0)
      xi = (/x1p,x2p/)
      xj = (/x1p,x2p/)
      fxinfo = (/fxc,fxcr,fxnp/)
      call rget(r,drdx,fxnvec,reps,xj,xi,fxinfo,
     &  depl,dtime,temp,props,nprops)
      fx = dsqrt(r(1)*r(1) + r(2)*r(2))
      fxd = drdx(1,1)*drdx(2,2)-drdx(2,1)*drdx(1,2)
      icount = 0
      newtmax=15
      if(xi(1).eq.(one))then
      newtmax = 50
      endif
      do while((icount.lt.newtmax).and.(dabs(fx).ge.toler))
      icount = icount+1
      if(dabs(fxd).gt.zero)then
      dxdr = reshape((/drdx(2,2),-drdx(2,1),
     &      -drdx(1,2),drdx(1,1)/),(/2,2/))/fxd
      xjupd = reshape(matmul(dxdr,reshape(r,(/2,1/))),(/2/))
      xj=xj-xjupd
      xj = (/max(dabs(xj(1)),x10),max(dabs(xj(2)),x20)/)
      fxinfo = (/fxc,fxcr,fxnp/)
      call rget(r,drdx,fxnvec,reps,xj,xi,fxinfo,
     &  depl,dtime,temp,props,nprops)
      fx = dsqrt(r(1)*r(1) + r(2)*r(2))
      fxd = drdx(1,1)*drdx(2,2)-drdx(2,1)*drdx(1,2)
      else
      xj = (/x1p,x2p/)
      fx = zero
      endif
      enddo
      x1 = max(xj(1),x10)
      x2 = max(xj(2),x20)
      fxn = min(dabs(fxnvec(1)),one)
c
c  add fxc contribution to Gamma 
      rx0 = one/(fxc-fxn)
      dxdfxc0 = dfxcde*(dtime*fxcr*rx0*rx0-rx0)
      dx1dfxc = dxdfxc0*x1
      dx2dfxc = dxdfxc0*x2      
      reps=reps+(/dx1dfxc,dx2dfxc/)
c     
      if(fxn.le.(1.d-3))then
      checkrx = .false.
      endif
      fxnr = fxnvec(2)
      dfxndx1 = fxnvec(3)
      dfxndx2 = fxnvec(4)
      dfxndfxc = fxnvec(5)
      
      xepl = xeplp*fxcp/fxc+depl
c
      tempstatev(lstskip+1) = xepl
      tempstatev(lstskip+2) = x1
      tempstatev(lstskip+3) = x2
      tempstatev(lstskip+4) = fxn
c
      x1eq = x1eq + x1*(fxc-fxn)
      x2eq = x2eq + x2*(fxc-fxn)
      plastic = plastic + xepl*(fxc-fxn)
c
      if(dabs(fxd).gt.0)then
      dx1de = dxdr(1,1)*reps(1)+dxdr(1,2)*reps(2)
      dx2de = dxdr(2,1)*reps(1)+dxdr(2,2)*reps(2)
      dfxnde = dfxndx1*dx1de+dfxndx2*dx2de+dfxndfxc*dfxcde
      
      
      dx1eqde = dx1eqde + dx1de*(fxc-fxn) +
     &      x1*(dfxcde - dfxnde)
      dfxcde = dfxnde
      fxc = fxn
      fxcp = fxnp
      fxcr = fxnr
      endif
c      endif
      ixvf = ixvf+1
      end do
c
      tempstatev(1) = plastic
      tempstatev(2) = x1eq
      tempstatev(3) = x2eq
c
      sqx1 = dsqrt(x1eq)
      sec = sig0*sqx1
      sy = siga + emusf*sfe*sec
c partial derivatives
c     d(sec)/d(epl)
      dsecdepl = half*sig0*dx1eqde/sqx1
c     d(sfe)/d(epl)
      dsfedepl = (sfe0*(one-sfel**(one/qe))**(one/pe-one)*
     &     sfel**(one/qe-one)/(pe*qe*rate))/dtime
c     total
      dsy = emusf*(sfe*dsecdepl+dsfedepl*sec)   
      return
      end
      