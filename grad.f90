program grad
	real*8 x(2),t,g(2),A(2,2),f,ans,aro(2)
	data x /2,2/
	data A /2,0,0,50/
	open(13,file="grad")
	write(13,"(5A13)")"f","x1","x2","g1","g2"
	do
	g=matmul(A,x)
	f=x(1)*x(1)+25.0*x(2)*x(2)
	write(13,"(5f13.6)")f,x,g
	call onedimenmul(g,g,2,ans)
	if(ans<0.01) exit
	t=ans
	aro=matmul(g,A)
	call onedimenmul(aro,g,2,ans)
	t=t/ans
	
	x=x-t*g
	end do
	close(13)
end program grad

	subroutine onedimenmul(m1,m2,n,ans)
	integer n
	real*8 m1(n),m2(n),ans
	ans=0
	do i=1,n
		ans=m1(i)*m2(i)+ans
	end do
	end subroutine onedimenmul
