
!1d lattice(2m+1) particles, all particles connected with each other in harmonic potential
!x(i)=displacement and momentum of ith particle in lattice starting from left end(i=1 for first particle in left size)


	module func
	implicit none
	contains


		double precision function f(k,a,q,y,z) !f=dp/dt=Force,force on nth particle={f(k,x[n-1],x[n],x[n+1])}
								        !dp/dt=-k((x(n)-x(n-1)+(x(n)-x(n+1)),   k=spring constant(analogous)
		double precision::q,y,z
		real::k,a
		f=-k*((y-q)+(y-z))+a*((z-y)**2-(y-q)**2)
		end function f



		double precision function g(q)       !g=dx/dt=p/mass(mass=1 here),p=momentum(hamiltons equation)
		double precision::q
		g=q
		end function g
	end module func

program pin
use func


implicit none

double precision::x(-4:1000),p(-4:1000),x0(-4:1000),x1(-4:1000),x2(-4:1000),x3(-4:1000),p0(-4:1000),p1(-4:1000),p2(-4:1000),&
p3(-4:1000),f0(-4:1000),f1(-4:1000),f2(-4:1000),f3(-4:1000)
real::t,tf,k,a,q,y,z,s,h,mass
integer::i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,m,i,i16,i17,ti,di,mode,pi,qi




mass=1.0
k=1.0
!stepsize=h
m=15



open(1,file="fer.dat")


do i16=-4,0
x0(i16)=0
x1(i16)=0
x2(i16)=0
x3(i16)=0
p0(i16)=0
p1(i16)=0
p2(i16)=0
p3(i16)=0
f0(i16)=0
f1(i16)=0
f2(i16)=0
f3(i16)=0
end do

do i16=2*m+1,2*m+5
x0(i16)=0
x1(i16)=0
x2(i16)=0
x3(i16)=0
p0(i16)=0
p1(i16)=0
p2(i16)=0
p3(i16)=0
f0(i16)=0
f1(i16)=0
f2(i16)=0
f3(i16)=0
end do




!we define m as the serial number of middlemost particle on lattice staring from left,(2m+1) particles as a whole 
!tf=time of observation



	write(*,*) "time of observation,nonlinearity parameter,initial mode of excitation"
	read(*,*) tf,a,mode



			

				do i=-4,2*m+5   !define initial momentum(zero for all), & endpoint particles cant move altogether as per question(fixed)
					p(i)=0
				end do





			
				do i1=-4,-1  !define initial displacement of particles(gaussian tilt in middle zone)

					x(i1)=0	
				end do	
	
				do i2=2*m+2,2*m+5 
					x(i2)=0	
				end do	
	
				do i3=0,2*m+1
					x(i3)=sin(real(mode)*3.1415926535*real(i3)/(2*m+1))  !1st mode

				end do	


	t=0
h=0.0001
	do while(t.le.tf)  !stop at time of observation




do i9=0,2*m+1
ti=t/h
pi=1.0/h
if(mod(ti,pi).eq.0) then
di=ti*h
	write(1,*) di,i9,x(i9)
end if
	end do




		do i4=1,2*m !start from 1, not zero, as end particles' position, momentum both are fixed at zero, so no need to change and update them with time goes on





!first runge kutta coefficient calculated for all ith particle's momentum change(p(i)->p(i)+(k1+2k2+2k3+k4)*h/6,ki depend supon force f(k,a,...)

!k1(-1),k1(2m+2)... are specified as they come up in the calculations for k4, they are essentially zero as end particles cant move.



			do i5=i4-4,i4+4    !started with i5=2, not 1, as, end particles are fixed, same for 2m+1th particle. k(1),k(2m+1) zero for same reason 
				
				if(i5.ge.1 .and. i5.le.2*m) then		
				x0(i5)=x(i5)
				x0(i5-1)=x(i5-1)
				x0(i5+1)=x(i5+1)
				p0(i5)=p(i5)	
				f0(i5)=f(k,a,x0(i5-1),x0(i5),x0(i5+1))
				else
				x0(i5)=0
				p0(i5)=0
				f0(i5)=0
				end if
			end do

		
			do i5=i4-4,i4+4    
				
				if(i5.ge.1 .and. i5.le.2*m) then
				x1(i5)=x0(i5)+h*p0(i5)/2
				x1(i5-1)=x0(i5-1)+h*p0(i5-1)/2
				x1(i5+1)=x0(i5+1)+h*p0(i5+1)/2
				p1(i5)=p0(i5)+h*f0(i5)/2	
				f1(i5)=f(k,a,x1(i5-1),x1(i5),x1(i5+1))
				else
				x1(i5)=0
				p1(i5)=0
				f1(i5)=0
				end if
											
			end do


			do i5=i4-4,i4+4    
				
				if(i5.ge.1 .and. i5.le.2*m) then
				x2(i5)=x0(i5)+h*p1(i5)/2
				x2(i5-1)=x0(i5-1)+h*p1(i5-1)/2
				x2(i5+1)=x0(i5+1)+h*p1(i5+1)/2
				p2(i5)=p0(i5)+h*f1(i5)/2	
				f2(i5)=f(k,a,x2(i5-1),x2(i5),x2(i5+1))
				else
				x2(i5)=0
				p2(i5)=0
				f2(i5)=0
				end if
											
			end do


			do i5=i4-4,i4+4    
				
				if(i5.ge.1 .and. i5.le.2*m) then
				x3(i5)=x0(i5)+h*p2(i5)
				x3(i5-1)=x0(i5-1)+h*p2(i5-1)
				x3(i5+1)=x0(i5+1)+h*p2(i5+1)
				p3(i5)=p0(i5)+h*f2(i5)	
				f3(i5)=f(k,a,x3(i5-1),x3(i5),x3(i5+1))
				else
				x3(i5)=0
				p3(i5)=0
				f3(i5)=0
				end if
											
			end do

		
			


		

!runge kutta coefficients for displacement change, dont need an array because depends upon one variable only(x(i)->x(i)+func(p(i))
			if(i4.ge.1 .and. i4.le.2*m) then			
			
			x(i4)=x(i4)+h*(p0(i4)+2*p1(i4)+2*p2(i4)+p3(i4))/(6*mass) 		!final runge kutta result(new x(i) and p(i))
			p(i4)=p(i4)+h*(f0(i4)+2*f1(i4)+2*f2(i4)+f3(i4))/6
			end if
		end do !for all particles(i4) at time t, again go back with new x(i) and p(i)s for t=t+2h and so on
t=t+h

	end do !for time=tf=time of observation



	!write to file

		

	

		close(1)

end program pin
		

	

















