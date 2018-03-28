!f在det=100时，迭代会变大，，，
!改变思路，为求f极小，试求它导数的平方的极小
!发现与原问题几乎等价，仍用原方法，想办法弄明白为什么DET＝100时沿G方向是递增的
program guangyichengzi
	implicit none
	integer i,j
	real*8 e,g2,f,g2ofg
	real*8 x(4),g(4),det(2),x1(4),gofg(4)!x(4) stand for x(1),x(2),x(3),x(4),
	real*8 t1,t2
	data det/100d0,100d0/
	data x /1.0d0,1d0,1d0,1d0/!data赋值好像不能在执行语句后(它与声明语句位置要求相同），尝试在循环语句中赋值失败
	data g /0d0,0d0,0d0,0d0/
	e=0.00001d0
	do
		t1=0d0;t2=0d0
		g2=0d0
		f=0d0
		do j=1,4
			x(j)=1d0
			g(j)=0d0
		end do
		
!		do i=1,2
!			call get_gofg(gofg,g2ofg,x,det)
!			if(g2ofg<=e)exit 
!			call jiabutansuo_ofg(x,gofg,det,t1,t2)
!			call duifenfa_ofg(t1,t2,x,det,gofg)
!		!call get_f(x,det,f)
!		print *,i
!		end do
	
!		if((((x(1)+x(2)-1)**2)<e.and.(x(2)>=1-e)).and.(g2<=e)) exit
!	!	if(((x(1)+x(2)-1)**2)>e) det(1)=10d0*det(1)
!	!	if(x(2)<=1d0-e) det(2)=10d0*det(2)
!		det=10*det
!	end do



		do i=1,2
			call get_g(g,g2,x,det)
			if(g2<=e)exit 
			call jiabutansuo(x,g,det,t1,t2)
			call duifenfa(t1,t2,x,det,g)
	!test
	t1=1D-7
	print *,t1
	call get_g(g,g2,x,det)
	print *,"x","f"
	do j=0,50
		x1=0.01d0*j*g+x
		call get_f(x1,det,f)
		print *,x1,f
	end do
	!test
		!call get_f(x,det,f)
		print *,i
		end do
	
		if((((x(1)+x(2)-1)**2)<e.and.(x(2)>=1-e)).and.(g2<=e)) exit
	!	if(((x(1)+x(2)-1)**2)>e) det(1)=10d0*det(1)
	!	if(x(2)<=1d0-e) det(2)=10d0*det(2)
		det=10*det
	end do
print *,x
end program guangyichengzi


subroutine get_g(g,g2,x,det)
	implicit none
	integer i,j
	real*8 x(4),det(2),g(4),g2
	real*8 u
	if(x(2)>=1d0)then
		u=0d0
	else
		u=1d0
	end if
	g(1)=-1d0+x(3)+2d0*det(1)*(x(1)+x(2)-1d0)
	g(2)=1d0+x(3)+2d0*det(1)*(x(1)+x(2)-1d0)+u*x(4)+2d0*det(2)*u*(x(2)-1d0)
	g(3)=x(1)+x(2)-1d0
	g(4)=u*(x(2)-1d0)
	g=-g
	g2=g(1)**2d0+g(2)**2d0+g(3)**2d0+g(4)**2d0
end subroutine 

subroutine get_f(x,det,f)
	implicit none
	real*8 x(4),det(2),f
	integer i,j
	real*8 u
	if(x(2)>=1d0)then
		u=0d0
	else
		u=(x(2)-1d0)
	end if
	f=(-x(1)+x(2))+x(3)*(x(1)+x(2)-1d0)+det(1)*(x(1)+x(2)-1d0)**2d0+x(4)*u+det(2)*u**2d0
end subroutine get_f

subroutine jiabutansuo(x,g,det,t1,t2)
	implicit none
	real*8 x(4),g(4),det(2),t1,t2
	integer i,j
	real*8 x1(4),f,f1
	t1=0d0
	t2=1d0
	call get_f(x,det,f)
!	i=1
!	do
!		t1=t1*i
!		x1=t1*g+x
!		call get_f(x1,det,f1)
!		if(f1>f)exit
!		i=i+1
!	end do
	i=1
	do
		t2=t2*i
		x1=t2*g+x
		call get_f(x1,det,f1)
		if(f1>f)exit
		i=i+1
	end do
end subroutine jiabutansuo

subroutine duifenfa(t1,t2,x,det,g)
	implicit none
	real*8 t1,t2,x(4),det(2),g(4)
	real*8, parameter::e=0.00001
	integer i,j
	real*8 x1(4),g1(4),f_t1,f_t2,f_t,g2 
	j=1
	do
		call get_ft(t1,x,det,f_t1)
		call get_ft(t2,x,det,f_t2)
		call get_ft((t1+t2)/2d0,x,det,f_t)
		if ((f_t)**2<e) exit
		if(f_t1*f_t<=0d0)then
			t2=(t1+t2)/2d0
		else
			t1=(t1+t2)/2d0
		end if
	j=j+1
	end do
	x=x+(t1+t2)/2d0*g
end subroutine

subroutine get_ft(t,x,det,f_t)
	implicit none
	integer i
	real*8 x(4),t,det(2),f_t
	real*8 g1(4),g2(4),g,dt,x1(4)
	call get_g(g1,g,x,det)
	call get_g(g2,g,g1*t+x,det)
	f_t=g1(1)*g2(1)+g1(2)*g2(2)+g1(3)*g2(3)+g1(4)*g2(4)
!	dt=1d-5
!	x1=x+(t-dt/2d0)*g
!	call get_f(x1,det,f1)
!	x1=x+(t+dt/2d0)*g
!	call get_f(x1,det,f2)
!	f_t=(f2-f1)/dt
end subroutine

subroutine  get_gofg(gofg,g2ofg,x,det)  
	implicit none
	integer i,j
	real*8 x(4),det(2),gofg(4),g2ofg
	real*8 u,g(4),g2
	if(x(2)>=1d0)then
		u=0d0
	else
		u=1d0
	end if
	call get_g(g,g2,x,det)
	g(1)=-1d0+x(3)+2d0*det(1)*(x(1)+x(2)-1d0)
	g(2)=1d0+x(3)+2d0*det(1)*(x(1)+x(2)-1d0)+u*x(4)+2d0*det(2)*u*(x(2)-1d0)
	g(3)=x(1)+x(2)-1d0
	g(4)=u*(x(2)-1d0)
	g2=g(1)**2d0+g(2)**2d0+g(3)**2d0+g(4)**2d0
	gofg(1)=2*g(1)*2d0*det(1)+2*g(2)*2d0*det(1)+2*g(3)
	gofg(2)=2*g(1)*2d0*det(1)+2*g(2)*(2d0*det(1)+2d0*det(2)*u)
	gofg(3)=2*g(1)+2*g(2)
	gofg(4)=2*g(2)*u
	gofg=-gofg
	g2ofg=gofg(1)**2+gofg(2)**2+gofg(3)**2+gofg(4)**2
end subroutine
