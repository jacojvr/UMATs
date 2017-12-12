      subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl,
     & ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, dtemp,
     & predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, nprops,
     & coords, drot, pnewdt, celent, df0, df1, noel, npt, layer,
     & kspt, kstep, kinc)
c
      implicit real*8(a-h,o-z)
c
      character*8 cmname
c
      dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens),
     & ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     & predef(1), dpred(1), props(nprops), coords(3), drot(3, 3),
     & df0(3, 3), df1(3, 3), tempstatev(nstatv), alpha(6), flow(6),
     &  xi(6), xiprev(6)
c
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     & enumax=.4999d0, newton=10, toler=1.0d-6)
c
c elastic properties
      call shearmod(eg,enu,temp,props,nprops)
      eg2=two*eg
      eg3=three*eg
      elam=(ebulk3-eg2)/three
c
c elastic stiffness
      do k1=1, ndi
        do k2=1, ndi
          ddsdde(k2, k1)=elam
        end do
        ddsdde(k1, k1)=eg2+elam
      end do
      do k1=ndi+1, ntens
        ddsdde(k1, k1)=eg
      end do
c
c recover and rotate shift tensor
      call rotsig(statev(1), drot, alpha, 1, ndi, nshr)
c
c previous effective stress
      do k1=1,ntens
        xiprev(k1)=stress(k1)-alpha(k1)
      enddo
c
c current predictor stress
      do k1=1, ntens
        do k2=1, ntens
          stress(k2)=stress(k2)+ddsdde(k2, k1)*dstran(k1)
        end do
      end do
c
c get trial effective stress:
      xixiprev = zero
      do k1=1,ntens
        xi(k1)=stress(k1)-alpha(k1)
        xixiprev = xi(k1)*xiprev(k1)
      enddo
c
c if sign change, then call load reversal subroutine
      if(xixiprev.lt.toler)then
        call cyclreversal(statev,nstatv)
      endif
c
c assign temporary state variables
      do k1=1,nstatv
        tempstatev(k1)=statev(k1)
      enddo
c
c calculate equivalent von mises stress
      smises=(xi(1)-xi(2))**2+(xi(2)-xi(3))**2+(xi(3)-xi(1))**2
      do k1=ndi+1,ntens
        smises=smises+six*xi(k1)**2
      end do
      smises=sqrt(smises/two)
c
c check yield surface
      call fcombined(sy0,dsy,sb0,dsb,zero,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
c
c determine if actively yielding
      if (smises.gt.(one+toler)*sy0) then
c
c actively yielding
c separate the hydrostatic from the deviatoric stress
c calculate the flow direction
      shydro=(stress(1)+stress(2)+stress(3))/three
      do k1=1,ndi
        flow(k1)=(stress(k1)-alpha(k1)-shydro)/smises
      end do
      do k1=ndi+1,ntens
      flow(k1)=(stress(k1)-alpha(k1))/smises
      end do
c
c solve for equivalent von mises stress
c and equivalent plastic strain increment using newton iteration
      sy=sy0
      sb=sb0
      deqpl=zero
      do kewton=1, newton
        rhs=smises-eg3*deqpl-sy-sb+sb0
        deqpl=dabs(deqpl+rhs/(eg3+dsy+dsb))
        call fcombined(sy,dsy,sb,dsb,deqpl,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
        if(abs(rhs).lt.toler*sy0) goto 10
      end do
c
c write warning message to .msg file
      write(7,2) newton
 2    format(//,30x,'***warning - plasticity algorithm did not ',
     & 'converge after ',i3,' iterations')
 10   continue
c
c update stress
      do k1=1,ndi
        alpha(k1)=alpha(k1)+(sb-sb0)*flow(k1)
        stress(k1)=alpha(k1)+flow(k1)*sy+shydro
      end do
      do k1=ndi+1,ntens
        alpha(k1)=alpha(k1)+(sb-sb0)*flow(k1)
        stress(k1)=alpha(k1)+flow(k1)*sy
      end do
c
c formulate the jacobian (material tangent)
c first calculate effective moduli
      effg=eg*(sy+sb-sb0)/smises
      effg2=two*effg
      effg3=three/two*effg2
      efflam=(ebulk3-effg2)/three
      effhrd=eg3*(dsy+dsb)/(eg3+dsy+dsb)-effg3
      do k1=1, ndi
        do k2=1, ndi
          ddsdde(k2, k1)=efflam
        end do
        ddsdde(k1, k1)=effg2+efflam
      end do
      do k1=ndi+1, ntens
        ddsdde(k1, k1)=effg
      end do
      do k1=1, ntens
        do k2=1, ntens
          ddsdde(k2, k1)=ddsdde(k2, k1)+effhrd*flow(k2)*flow(k1)
        end do
      end do
      endif
c
c update state variable array
      do k1=1,ntens
        tempstatev(k1)=alpha(k1)
      enddo
      do k1=1,nstatv
        statev(k1)=tempstatev(k1)
      enddo
c
      return
      end