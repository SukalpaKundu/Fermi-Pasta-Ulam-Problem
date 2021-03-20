






implicit none


double precision::p(200),q(200),p0(200),q0(200),f0(200),p1(200),q1(200),f1(200),p2(200),q2(200),f2(200),&
p3(200),q3(200),f3(200),h,tf,s,e,a1,b1,t,a(200,200,200),energy
integer::n,i,mode,l,m,j,ti,di,pi
open(1,file='a.dat')



write(*,*) "tf,initial energy,,nonlinearity(e),mode"
read(*,*) tf,energy,e,mode


n=32
h=0.1



!initial normal mode description

do i=1,n
q(i)=0
end do



q(mode)=sqrt(2*energy)/(2*sin(dble(mode)*3.1419/(2*(dble(n)+1))))




!initial momentum to each normal mode

do i=1,n
p(i)=0
end do






do i=1,n
	do l=1,n
		do m=1,n

		if(i.eq.l+m .or. i.eq.-l+m .or. i.eq.l-m) then
		a1=1
		else 
		a1=0
		end if



		if(i+l+m.eq.n+n+2 .or. i-l+m.eq.n+n+2 .or. i+l-m.eq.n+n+2 .or. i-l-m.eq.n+n+2) then
		b1=1
		else
		b1=0
		end if
	
	a(i,l,m)=(a1-b1)*2*sin(dble(i)*3.1419/(2*(dble(n)+1)))*2*sin(dble(l)&
*3.1419/(2*(dble(n)+1)))*2*sin(dble(m)*3.1419/(2*(dble(n)+1)))/sqrt(2*(dble(n)+1))

end do
	end do
		end do











t=1
do while(t .le. tf)



		!runge kutta coefficients
		
			do i=1,n
			do j=1,n
					
				q0(j)=q(j)
				p0(j)=p(j)
			end do
				s=0
				forall(l=1:n,m=1:n)

					s=s-a(i,l,m)*q0(l)*q0(m)*e

			end forall
				
				
				
				f0(i)=s-(2*sin(dble(i)*3.1419/(2*(dble(n)+1))))**2*q0(i)
		end do
		


		
		
			do i=1,n
			do j=1,n
				
				q1(j)=q0(j)+h*p0(j)/2
				p1(j)=p0(j)+h*f0(j)/2
			end do
				s=0
				forall(l=1:n,m=1:n)

					s=s-a(i,l,m)*q1(l)*q1(m)*e
					end forall
				
				f1(i)=s-(2*sin(dble(i)*3.1419/(2*(dble(n)+1))))**2*q1(i)
		end do
		

		
		do i=1,n
		
			do j=1,n
					
				q2(j)=q0(j)+h*p1(j)/2
				p2(j)=p0(j)+h*f1(j)/2
			end do
				s=0
				forall(l=1:n,m=1:n)
					s=s-a(i,l,m)*q2(l)*q2(m)*e

					end forall
				f2(i)=s-(2*sin(dble(i)*3.1419/(2*(dble(n)+1))))**2*q2(i)
		end do


		
			
		do i=1,n
			do j=1,n
				
				q3(j)=q0(j)+h*p2(j)
				p3(j)=p0(j)+h*f2(j)
			end do
				s=0
				forall(l=1:n,m=1:n)


					s=s-a(i,l,m)*q3(l)*q3(m)*e

					end forall
				
				f3(i)=s-(2*sin(dble(i)*3.1419/(2*(dble(n)+1))))**2*q3(i)
		end do
	 !ends calculation for rk coefficients




!new ith variable
do i=1,n
q(i)=q(i)+h*(p0(i)+2*p1(i)+2*p2(i)+p3(i))/6
p(i)=p(i)+h*(f0(i)+2*f1(i)+2*f2(i)+f3(i))/6


end do

t=t+h	


		!write

	ti=t/h
	pi=10/h
	if(mod(ti,pi).eq.0 .and. t.ge.1) then

			
		write(1,*) t,0.5*((2*sin(dble(mode)*3.1419/(2*(dble(n)+1))))**2*abs(q(mode))**2+abs(p(mode))**2)
		
	end if

ti=t/h
	pi=100/h
	if(mod(ti,pi).eq.0 .and. t.ge.1) then
	
		
		write(*,*) t
	end if


end do !for time
close(1)

end












			
























