      subroutine rget(r,drdxj,fxnvec,reps,xj,xi,fxinfo,
     &  depl,dtime,temp,props,nprops)
c
      implicit real*8(a-h,o-z)
      dimension props(nprops),r(2),drdxj(2,2),fxinfo(3),reps(2),
     &  xj(2),xi(2),fxnvec(5)
c
      parameter(zero=0.d0,half=0.5d0,one=1.d0,two=2.d0,
     & toler=1.d-10,ratelim=1.d-8)
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
c evolution of misorient:      
      cld = props(11)
c stage 4
      cg = props(12)
      rg = props(13)
c storage
      c1 = props(14)
c dynamic recovery:
      c20 = props(15)
      a02 = props(16)
      rate02 = props(17)
c thermal recovery
      c30 = props(18)
      r3 = props(19)
      a03 = props(20)
c recrystallisation
      cx0 = props(21)
      a0x = props(22)
      cxl = props(23)
      rxl = props(24)
      rxa = props(25)
      rxb = props(26)
      cxc = props(27)
c
c get info from previous values: contained in fxinfo:
      fxc = fxinfo(1)
      fxcr = dabs(fxinfo(2))
      fxnp = fxinfo(3)
c
      rate = depl/dtime
      if(rate.lt.ratelim)then
      rate = ratelim
      endif
c
      if(temp.gt.et0)then
      emu = emu0 - ed0/(dexp(et0/temp)-one)
      r2m = -temp/(emu*a02)
      r3c = c30*dexp(-a03/temp)*dtime
      r5c0 = cx0*dexp(-a0x/temp)*emu*dtime
      else
c use constants:
      emu = emu0
      r2m = -a02
      r3c = c30*dtime
      r5c0 = cx0*emu*dtime
      endif
c isv's at previous convergence and current guess:
      x1p = max(xi(1),one)
      x2p = max(xi(2),zero)
      x1 = max(xj(1),one)
      x2 = max(xj(2),zero)
c
c growth of next recrystallised volume:
      cldbar = min(x2,one)
      r5c1 = (one-dexp(-cxl*(cldbar)**rxl))
      r5c = r5c0*r5c1*(x1)
c interfacial area 
      fxn = fxnp
      rxg = fxc*((fxn/fxc)**rxa)*((one-fxn/fxc)**rxb)*
     & (one+cxc*(one-fxc))
      drxg = rxa*((fxn/fxc)**(rxa-one))*
     & ((one-fxn/fxc)**rxb)*(one+cxc*(one-fxc)) -
     & rxb*((fxn/fxc)**rxa)*((one-fxn/fxc)**(rxb-one))*
     & (one+cxc*(one-fxc))
      fxnr = dabs(r5c*rxg)
      ffxn = fxn - fxnp - fxnr
c resolve residual
      icount = 0
      do while((icount.lt.15).and.(dabs(ffxn).gt.toler))
      icount = icount+1
      dffxn = one - half*r5c*drxg
      if(dabs(dffxn).lt.toler)then
      dffxn = toler
      endif
      fxn = min(dabs(fxn - ffxn/dffxn),fxc-toler)
      rxg = fxc*((fxn/fxc)**rxa)*((one-fxn/fxc)**rxb)*
     & (one+cxc*(one-fxc))
      drxg = rxa*((fxn/fxc)**(rxa-one))*
     & ((one-fxn/fxc)**rxb)*(one+cxc*(one-fxc)) -
     & rxb*((fxn/fxc)**rxa)*((one-fxn/fxc)**(rxb-one))*
     & (one+cxc*(one-fxc))
      fxnr = dabs(r5c*rxg)
      ffxn = fxn - fxnp - fxnr
      end do
c partial : change of fxn with respect to x1 and x2:
c partial gradients d(fxn)/d(x1)   
      ddfxdr = (one-r5c*drxg)
      if(dabs(ddfxdr).lt.toler)then
      ddfxdr = toler
      endif 
      dfxndx1 = (r5c0*r5c1*rxg)/ddfxdr
      dmdx2 = zero
      if(x2.lt.one)then
      dmdx2 = one
      endif
      dr5c1dm = rxl*cxl*dexp(-cxl*(cldbar)**rxl)*
     & (cldbar)**(rxl-one)
      dfxndm = (r5c0*dr5c1dm*x1*rxg)/ddfxdr
c partial gradients d(fxn)/d(x2)
      dfxndx2 = dfxndm*dmdx2
c partial d(fxn)/d(fxc)
      drxgdfxc = ((fxn/fxc)**rxa)*((one-fxn/fxc)**rxb)*
     & (one+cxc*(one-fxc)) 
     & - rxa*((fxn/fxc)**rxa)*((one-fxn/fxc)**rxb)*
     & (one+cxc*(one-fxc))
     & + rxb*((fxn/fxc)**(rxa+one))*((one-fxn/fxc)**(rxb-one))*
     & (one+cxc*(one-fxc))
     & - cxc*fxc*((fxn/fxc)**rxa)*((one-fxn/fxc)**rxb)
      dfxndfxc = r5c*drxgdfxc
c residual equations on the isv values:
      c2 = c20*(rate/rate02)**r2m
      dc2de = c20*r2m*((rate/rate02)**(r2m-one))/
     & (dtime*rate02)
      sqx1p = dsqrt(x1p)
      sqx1 = dsqrt(x1)
      hx2 = cld
      hx1temp = -r3c*(x1**r3+x1p**r3)
      hx1 = (cg)*(x2)**rg+c1*sqx1-c2*x1
      if(fxn.ge.fxc)then
      rx0 = zero
      else
      rx0 = one/(fxc - fxn)
      endif
      rfxc = fxcr*rx0
      drx = rfxc*rx0
      r1 = x1-x1p-hx1*depl-hx1temp+x1*rfxc
      r2 = x2-x2p-hx2*depl+x2*rfxc
      r = (/ r1 , r2 /)
c partial : change in residual with respect to x1 and x2:
      dhrdrtemp=-r3*r3c*x1**(r3-one)
      dhrdr=half*c1/sqx1-c2
c partial gradients d(fr1)/d(x1)
      dr1dr0=one-dhrdr*depl-dhrdrtemp+rfxc
      dr1dx1=dr1dr0+x1*drx*dfxndx1
c partial gradients d(fr1)/d(x2)
      dr1dx2=-depl*rg*(cg)*(x2)**(rg-one)+
     & x1*drx*dfxndx2
c partial gradients d(f2)/d(x1)
      dr2dx1 = x2*drx*dfxndx1
c partial gradients d(f2)/d(x2)
      dr2dl0= one + rfxc
      dr2dx2 = dr2dl0 + x2*drx*dfxndx2
      drdxj=reshape((/dr1dx1,dr2dx1,dr1dx2,dr2dx2/),(/2,2/))
c
      dr1de = hx1 - depl*dc2de*x1
      dr2de = hx2
      reps = (/dr1de,dr2de/)  
c exchange new supplementary info using fxinfo:
      fxnvec(1) = fxn
      fxnvec(2) = fxnr
      fxnvec(3) = dfxndx1
      fxnvec(4) = dfxndx2
      fxnvec(5) = dfxndfxc
c      
      return
      end