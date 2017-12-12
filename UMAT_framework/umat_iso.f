      subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl,
     & ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, dtemp,
     & predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, nprops,
     & coords, drot, pnewdt, celent, df0, df1, noel, npt, layer,
     & kspt, kstep, kinc)
c
      implicit real*8(a-h,o-z)
      character*8 cmname
c
      dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens),
     & ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     & predef(1), dpred(1), props(nprops), coords(3), drot(3, 3),
     & df0(3, 3), df1(3, 3), flow(6), tempstatv(nstatv)
c
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     & enumax=.4999d0, newton=10, toler=1.0d-6)
c
c assign temporary state variables
      do k1=1,nstatv
       tempstatv(k1)=statev(k1)
      enddo
c
c elastic properties
      call shearmod(eg,enu,temp,props,nprops)
      eg2=two*eg
      eg3=three*eg
      emod = eg2*(1.d0+enu)
      ebulk3=emod/(one-two*enu)
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
c calculate predictor stress and elastic strain
      do k1=1, ntens
        do k2=1, ntens
          stress(k2)=stress(k2)+ddsdde(k2, k1)*dstran(k1)
        end do
      end do
c
c calculate equivalent von mises stress
      smises=(stress(1)-stress(2))**2+(stress(2)-stress(3))**2
     & +(stress(3)-stress(1))**2
      do k1=ndi+1,ntens
        smises=smises+six*stress(k1)**2
      end do
      smises=sqrt(smises/two)
c
      call fisotropic(sy,dsy,zero,dtime,temp,
     & statev,tempstatv,nstatv,props,nprops)
c
c determine if actively yielding
      if (smises.gt.(one+toler)*sy) then
c
c actively yielding
c separate the hydrostatic from the deviatoric stress
c calculate the flow direction
      shydro=(stress(1)+stress(2)+stress(3))/three
      do k1=1,ndi
        flow(k1)=(stress(k1)-shydro)/smises
      end do
      do k1=ndi+1, ntens
        flow(k1)=stress(k1)/smises
      end do
c
c solve for equivalent von mises stress
c and equivalent plastic strain increment using newton iteration
      deqpl=zero
      do kewton=1, newton
        rhs=smises-eg3*deqpl-sy
        deqpl=deqpl+rhs/(eg3+dsy)
	call fisotropic(sy,dsy,deqpl,dtime,temp,
     & statev,tempstatv,nstatv,props,nprops)
        if(abs(rhs).lt.toler) goto 10
      end do
c
c write warning message to .msg file
      write(7,2) newton
 2    format(//,30x,'***warning - plasticity algorithm did not ',
     & 'converge after ',i3,' iterations')
 10   continue
c
c update stress, elastic and plastic strains and
c equivalent plastic strain
      do k1=1,ndi
        stress(k1)=flow(k1)*sy+shydro
      end do
      do k1=ndi+1,ntens
        stress(k1)=flow(k1)*sy
      end do
c
c formulate the jacobian (material tangent)
c first calculate effective moduli
      effg=eg*sy/smises
      effg2=two*effg
      effg3=three/two*effg2
      efflam=(ebulk3-effg2)/three
      effhrd=eg3*dsy/(eg3+dsy)-effg3
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
      do k1=1,nstatv
       statev(k1)=tempstatv(k1)
      enddo
c
      return
      end
c