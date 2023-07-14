#####################################################
#####################################################
using Plots
using Interpolations
using QuadGK
using Printf
#####################################################
#####################################################
function getpaths(n1,n2,n3,n4) #generate four paths that define the c domain
#paths defined in counterclockwise orientation for Green's Thm

paths=Dict();

t=collect(range(pi/4,2*pi,n1)); #outer circle of radius 3
x=3*cos.(t)
y=3*sin.(t)
p1=Matrix{Float64}(undef,n1,2)
p1[:,1]=x
p1[:,2]=y
paths["p1"]=p1;

x=collect(range(3,1,n2));  #horizontal line
p2=Matrix{Float64}(undef,n2,2)
p2[:,1]=x
p2[:,2].=0.0
paths["p2"]=p2;

t=t[end:-1:1] #inner circle of radius one
t=collect(range(2*pi,pi/4,n3));
x=cos.(t)
y=sin.(t)
p3=Matrix{Float64}(undef,n3,2)
p3[:,1]=x
p3[:,2]=y
paths["p3"]=p3;

x=collect(range(cos(pi/4),3*cos(pi/4),n4)); #line of slope 1
p4=Matrix{Float64}(undef,n4,2)
p4[:,1]=x
p4[:,2]=x
paths["p4"]=p4;

return paths

end
################################################################
################################################################

function greens7(nn,n)
#nn points for each path
#n segments in x and y for ODE method to get integrand
#Actual integrand for this example is f(x,y)=x^2+y^2
n1=nn
n2=nn
n3=nn
n4=nn
nx = n
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for ODE method
z = Array{Float64}(undef, nx, ny)

#rk4 method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  k1=f(x1[i-1],y1[j])
  k2=f((x1[i-1]+x1[i])/2,y1[j])
  k4=f(x1[i],y1[j])
  z[j,i]=z[j,i-1]+(hx/6)*(k1+4*k2+k4)
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

#interpolate outer circle
t1 = 0:h1:1
###approximate dy at knots
dy=Array{Float64}(undef,n1)
dy[1]=(-(25/12)*p1[1,2]+4*p1[2,2]-3*p1[3,2]+(4/3)*p1[4,2]-(1/4)*p1[5,2])/h1
dy[2]=(-(25/12)*p1[2,2]+4*p1[3,2]-3*p1[4,2]+(4/3)*p1[5,2]-(1/4)*p1[6,2])/h1
for i=3:n1-2
dy[i]=(-p1[i+2,2]+8*p1[i+1,2]-8*p1[i-1,2]+p1[i-2,2])/(12*h1)
end
dy[n1-1]=((25/12)*p1[n1-1,2]-4*p1[n1-2,2]+3*p1[n1-3,2]-(4/3)*p1[n1-4,2]+(1/4)*p1[n1-5,2])/h1
dy[n1]=((25/12)*p1[n1,2]-4*p1[n1-1,2]+3*p1[n1-2,2]-(4/3)*p1[n1-3,2]+(1/4)*p1[n1-4,2])/h1
int1=itp(p1[1,1],p1[1,2]).*dy[1]+itp(p1[n1,1],p1[n1,2]).*dy[n1]
for i=2:n1-1
int1+=2*itp(p1[i,1],p1[i,2]).*dy[i]
end
int1*=(h1/2)

#line integral for horizontal line
t2 = 0:h2:1
#approximate dy at knots
dy=Array{Float64}(undef,n2)
dy[1]=(-(25/12)*p2[1,2]+4*p2[2,2]-3*p2[3,2]+(4/3)*p2[4,2]-(1/4)*p2[5,2])/h2
dy[2]=(-(25/12)*p2[2,2]+4*p2[3,2]-3*p2[4,2]+(4/3)*p2[5,2]-(1/4)*p2[6,2])/h2
for i=3:n2-2
dy[i]=(-p2[i+2,2]+8*p2[i+1,2]-8*p2[i-1,2]+p2[i-2,2])/(12*h2)
end
dy[n2-1]=((25/12)*p2[n2-1,2]-4*p2[n2-2,2]+3*p2[n2-3,2]-(4/3)*p2[n2-4,2]+(1/4)*p2[n2-5,2])/h2
dy[n2]=((25/12)*p2[n2,2]-4*p2[n2-1,2]+3*p2[n2-2,2]-(4/3)*p2[n2-3,2]+(1/4)*p2[n2-4,2])/h2
int2=itp(p2[1,1],p2[1,2]).*dy[1]+itp(p2[n2,1],p2[n2,2]).*dy[n2]
for i=2:n2-1
int2+=2*itp(p2[i,1],p2[i,2]).*dy[i]
end
int2*=(h2/2)

#line integral for inner circle
t3 = 0:h3:1
#approximate dy at knots
dy=Array{Float64}(undef,n3)
dy[1]=(-(25/12)*p3[1,2]+4*p3[2,2]-3*p3[3,2]+(4/3)*p3[4,2]-(1/4)*p3[5,2])/h3
dy[2]=(-(25/12)*p3[2,2]+4*p3[3,2]-3*p3[4,2]+(4/3)*p3[5,2]-(1/4)*p3[6,2])/h3
for i=3:n3-2
dy[i]=(-p3[i+2,2]+8*p3[i+1,2]-8*p3[i-1,2]+p3[i-2,2])/(12*h3)
end
dy[n3-1]=((25/12)*p3[n3-1,2]-4*p3[n3-2,2]+3*p3[n3-3,2]-(4/3)*p3[n3-4,2]+(1/4)*p3[n3-5,2])/h3
dy[n3]=((25/12)*p3[n3,2]-4*p3[n3-1,2]+3*p3[n3-2,2]-(4/3)*p3[n3-3,2]+(1/4)*p3[n3-4,2])/h3
int3=itp(p3[1,1],p3[1,2]).*dy[1]+itp(p3[n3,1],p3[n3,2]).*dy[n3]
for i=2:n1-1
int3+=2*itp(p3[i,1],p3[i,2]).*dy[i]
end
int3*=(h3/2)



#line integral for line w/ slope 1
t4 = 0:h4:1
#approximate dy at knots
dy=Array{Float64}(undef,n4)
dy[1]=(-(25/12)*p4[1,2]+4*p4[2,2]-3*p4[3,2]+(4/3)*p4[4,2]-(1/4)*p4[5,2])/h4
dy[2]=(-(25/12)*p4[2,2]+4*p4[3,2]-3*p4[4,2]+(4/3)*p4[5,2]-(1/4)*p4[6,2])/h4
for i=3:n4-2
dy[i]=(-p4[i+2,2]+8*p4[i+1,2]-8*p4[i-1,2]+p4[i-2,2])/(12*h4)
end
dy[n4-1]=((25/12)*p4[n4-1,2]-4*p4[n4-2,2]+3*p4[n4-3,2]-(4/3)*p4[n4-4,2]+(1/4)*p4[n4-5,2])/h4
dy[n4]=((25/12)*p4[n4,2]-4*p4[n4-1,2]+3*p4[n4-2,2]-(4/3)*p4[n4-3,2]+(1/4)*p4[n4-4,2])/h4
int4=itp(p4[1,1],p4[1,2]).*dy[1]+itp(p4[n4,1],p4[n4,2]).*dy[n4]
for i=2:n4-1
int4+=2*itp(p4[i,1],p4[i,2]).*dy[i]
end
int4*=(h4/2)

int=(int1+int2+int3+int4)


exactarea=35*pi
err=abs((exactarea-int)/exactarea)
return err
end


################################################################
################################################################


################################################################
################################################################

function greens6(nn,n)
#nn points for each path
#n segments in x and y for ODE method to get integrand
#Actual integrand for this example is f(x,y)=x^2+y^2
n1=nn
n2=nn
n3=nn
n4=nn
nx = n
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for ODE method
z = Array{Float64}(undef, nx, ny)

#rk4 method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  k1=f(x1[i-1],y1[j])
  k2=f((x1[i-1]+x1[i])/2,y1[j])
  k4=f(x1[i],y1[j])
  z[j,i]=z[j,i-1]+(hx/6)*(k1+4*k2+k4)
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

#interpolate outer circle
t1 = 0:h1:1
itp1 = Interpolations.scale(interpolate(p1, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t1, 1:2)
###approximate dy at knots
dy=Array{Float64}(undef,n1)
dy[1]=(-(25/12)*p1[1,2]+4*p1[2,2]-3*p1[3,2]+(4/3)*p1[4,2]-(1/4)*p1[5,2])/h1
dy[2]=(-(25/12)*p1[2,2]+4*p1[3,2]-3*p1[4,2]+(4/3)*p1[5,2]-(1/4)*p1[6,2])/h1
for i=3:n1-2
dy[i]=(-p1[i+2,2]+8*p1[i+1,2]-8*p1[i-1,2]+p1[i-2,2])/(12*h1)
end
dy[n1-1]=((25/12)*p1[n1-1,2]-4*p1[n1-2,2]+3*p1[n1-3,2]-(4/3)*p1[n1-4,2]+(1/4)*p1[n1-5,2])/h1
dy[n1]=((25/12)*p1[n1,2]-4*p1[n1-1,2]+3*p1[n1-2,2]-(4/3)*p1[n1-3,2]+(1/4)*p1[n1-4,2])/h1
###interpolate dy
itpdy1 = linear_interpolation(t1, dy)
###approximate interal
int1, error = quadgk(t -> itp(itp1(t,1),itp1(t,2))*itpdy1(t), 0, 1)

#line integral for horizontal line
t2 = 0:h2:1
itp2 = Interpolations.scale(interpolate(p2, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t2, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n2)
dy[1]=(-(25/12)*p2[1,2]+4*p2[2,2]-3*p2[3,2]+(4/3)*p2[4,2]-(1/4)*p2[5,2])/h2
dy[2]=(-(25/12)*p2[2,2]+4*p2[3,2]-3*p2[4,2]+(4/3)*p2[5,2]-(1/4)*p2[6,2])/h2
for i=3:n2-2
dy[i]=(-p2[i+2,2]+8*p2[i+1,2]-8*p2[i-1,2]+p2[i-2,2])/(12*h2)
end
dy[n2-1]=((25/12)*p2[n2-1,2]-4*p2[n2-2,2]+3*p2[n2-3,2]-(4/3)*p2[n2-4,2]+(1/4)*p2[n2-5,2])/h2
dy[n2]=((25/12)*p2[n2,2]-4*p2[n2-1,2]+3*p2[n2-2,2]-(4/3)*p2[n2-3,2]+(1/4)*p2[n2-4,2])/h2
#interpolate dy
itpdy2 = linear_interpolation(t2, dy)
#approximate interal
int2, error = quadgk(t -> itp(itp2(t,1),itp2(t,2))*itpdy2(t), 0, 1)

#line integral for inner circle
t3 = 0:h3:1
itp3 = Interpolations.scale(interpolate(p3, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t3, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n3)
dy[1]=(-(25/12)*p3[1,2]+4*p3[2,2]-3*p3[3,2]+(4/3)*p3[4,2]-(1/4)*p3[5,2])/h3
dy[2]=(-(25/12)*p3[2,2]+4*p3[3,2]-3*p3[4,2]+(4/3)*p3[5,2]-(1/4)*p3[6,2])/h3
for i=3:n3-2
dy[i]=(-p3[i+2,2]+8*p3[i+1,2]-8*p3[i-1,2]+p3[i-2,2])/(12*h3)
end
dy[n3-1]=((25/12)*p3[n3-1,2]-4*p3[n3-2,2]+3*p3[n3-3,2]-(4/3)*p3[n3-4,2]+(1/4)*p3[n3-5,2])/h3
dy[n3]=((25/12)*p3[n3,2]-4*p3[n3-1,2]+3*p3[n3-2,2]-(4/3)*p3[n3-3,2]+(1/4)*p3[n3-4,2])/h3
#interpolate dy
itpdy3 = linear_interpolation(t3, dy)
#approximate interal
int3, error = quadgk(t -> itp(itp3(t,1),itp3(t,2))*itpdy3(t), 0, 1)



#line integral for line w/ slope 1
t4 = 0:h4:1
itp4 = Interpolations.scale(interpolate(p4, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t4, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n4)
dy[1]=(-(25/12)*p4[1,2]+4*p4[2,2]-3*p4[3,2]+(4/3)*p4[4,2]-(1/4)*p4[5,2])/h4
dy[2]=(-(25/12)*p4[2,2]+4*p4[3,2]-3*p4[4,2]+(4/3)*p4[5,2]-(1/4)*p4[6,2])/h4
for i=3:n4-2
dy[i]=(-p4[i+2,2]+8*p4[i+1,2]-8*p4[i-1,2]+p4[i-2,2])/(12*h4)
end
dy[n4-1]=((25/12)*p4[n4-1,2]-4*p4[n4-2,2]+3*p4[n4-3,2]-(4/3)*p4[n4-4,2]+(1/4)*p4[n4-5,2])/h4
dy[n4]=((25/12)*p4[n4,2]-4*p4[n4-1,2]+3*p4[n4-2,2]-(4/3)*p4[n4-3,2]+(1/4)*p4[n4-4,2])/h4
#interpolate dy
itpdy4 = linear_interpolation(t4, dy)
#approximate interal
int4, error = quadgk(t -> itp(itp4(t,1),itp4(t,2))*itpdy4(t), 0, 1)

int=(int1+int2+int3+int4)


exactarea=35*pi
err=abs((exactarea-int)/exactarea)
return err
end


################################################################
################################################################

function greens5(nn,n)
#nn points for each path
#n segments in x and y for ODE method to get integrand
#Actual integrand for this example is f(x,y)=x^2+y^2
n1=nn
n2=nn
n3=nn
n4=nn
#set parameters for ODE method
nx = n
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for ODE method
z = Array{Float64}(undef, nx, ny)

#rk4 method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  k1=f(x1[i-1],y1[j])
  k2=f((x1[i-1]+x1[i])/2,y1[j])
  k4=f(x1[i],y1[j])
  z[j,i]=z[j,i-1]+(hx/6)*(k1+4*k2+k4)
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

#interpolate outer circle
t1 = 0:h1:1
itp1 = Interpolations.scale(interpolate(p1, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t1, 1:2)
###approximate dy at knots
dy=Array{Float64}(undef,n1)
dy[1]=(-3*p1[1,2]+4*p1[2,2]-p1[3,2])/(2*h1)
for i=2:n1-1
dy[i]=(p1[i+1,2]-p1[i-1,2])/(2*h1)
end
dy[n1]=(3*p1[n1,2]-4*p1[n1-1,2]+p1[n1-2,2])/(2*h1)
###interpolate dy
itpdy1 = linear_interpolation(t1, dy)
###approximate interal
int1, error = quadgk(t -> itp(itp1(t,1),itp1(t,2))*itpdy1(t), 0, 1)

#line integral for horizontal line
t2 = 0:h2:1
itp2 = Interpolations.scale(interpolate(p2, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t2, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n2)
dy[1]=(-3*p2[1,2]+4*p2[2,2]-p2[3,2])/(2*h2)
for i=2:n2-1
dy[i]=(p2[i+1,2]-p2[i-1,2])/(2*h2)
end
dy[n2]=(3*p2[n2,2]-4*p2[n2-1,2]+p2[n2-2,2])/(2*h2)
#interpolate dy
itpdy2 = linear_interpolation(t2, dy)
#approximate interal
int2, error = quadgk(t -> itp(itp2(t,1),itp2(t,2))*itpdy2(t), 0, 1)

#line integral for inner circle
t3 = 0:h3:1
itp3 = Interpolations.scale(interpolate(p3, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t3, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n3)
dy[1]=(-3*p3[1,2]+4*p3[2,2]-p3[3,2])/(2*h3)
for i=2:n3-1
dy[i]=(p3[i+1,2]-p3[i-1,2])/(2*h3)
end
dy[n3]=(3*p3[n3,2]-4*p3[n3-1,2]+p3[n3-2,2])/(2*h3)
#interpolate dy
itpdy3 = linear_interpolation(t3, dy)
#approximate interal
int3, error = quadgk(t -> itp(itp3(t,1),itp3(t,2))*itpdy3(t), 0, 1)



#line integral for line w/ slope 1
t4 = 0:h4:1
itp4 = Interpolations.scale(interpolate(p4, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t4, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n4)
dy[1]=(-3*p4[1,2]+4*p4[2,2]-p4[3,2])/(2*h4)
for i=2:n4-1
dy[i]=(p4[i+1,2]-p4[i-1,2])/(2*h4)
end
dy[n4]=(3*p4[n4,2]-4*p4[n4-1,2]+p4[n4-2,2])/(2*h4)
#interpolate dy
itpdy4 = linear_interpolation(t4, dy)
#approximate interal
int4, error = quadgk(t -> itp(itp4(t,1),itp4(t,2))*itpdy4(t), 0, 1)

int=(int1+int2+int3+int4)

exactarea=35*pi
err=abs((exactarea-int)/exactarea)
return err
end
################################################################
################################################################

function greens4(nn,n)
#nn points for each path
#n segments in x and y for ODE method to get integrand
#Actual integrand for this example is f(x,y)=x^2+y^2
n1=nn
n2=nn
n3=nn
n4=nn
#set parameters for ODE method
nx = n
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for ODE method
z = Array{Float64}(undef, nx, ny)

#rk4 method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  k1=f(x1[i-1],y1[j])
  k2=f((x1[i-1]+x1[i])/2,y1[j])
  k4=f(x1[i],y1[j])
  z[j,i]=z[j,i-1]+(hx/6)*(k1+4*k2+k4)
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

#interpolate outer circle
t1 = 0:h1:1
itp1 = Interpolations.scale(interpolate(p1, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t1, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n1)
dy[1]=(p1[2,2]-p1[1,2])/h1
for i=2:n1
dy[i]=(p1[i,2]-p1[i-1,2])/h1
end
#interpolate dy
itpdy1 = linear_interpolation(t1, dy)
#approximate interal
int1, error = quadgk(t -> itp(itp1(t,1),itp1(t,2))*itpdy1(t), 0, 1)

#line integral for horizontal line
t2 = 0:h2:1
itp2 = Interpolations.scale(interpolate(p2, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t2, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n2)
dy[1]=(p2[2,2]-p2[1,2])/h2
for i=2:n2
dy[i]=(p2[i,2]-p2[i-1,2])/h2
end
#interpolate dy
itpdy2 = linear_interpolation(t2, dy)
#approximate interal
int2, error = quadgk(t -> itp(itp2(t,1),itp2(t,2))*itpdy2(t), 0, 1)

#line integral for inner circle
t3 = 0:h3:1
itp3 = Interpolations.scale(interpolate(p3, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t3, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n3)
dy[1]=(p3[2,2]-p3[1,2])/h3
for i=2:n3
dy[i]=(p3[i,2]-p3[i-1,2])/h3
end
#interpolate dy
itpdy3 = linear_interpolation(t3, dy)
#approximate interal
int3, error = quadgk(t -> itp(itp3(t,1),itp3(t,2))*itpdy3(t), 0, 1)



#line integral for line w/ slope 1
t4 = 0:h4:1
itp4 = Interpolations.scale(interpolate(p4, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t4, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n4)
dy[1]=(p4[2,2]-p4[1,2])/h4
for i=2:n4
dy[i]=(p4[i,2]-p4[i-1,2])/h4
end
#interpolate dy
itpdy4 = linear_interpolation(t4, dy)
#approximate interal
int4, error = quadgk(t -> itp(itp4(t,1),itp4(t,2))*itpdy4(t), 0, 1)

int=(int1+int2+int3+int4)

exactarea=35*pi
err=abs((exactarea-int)/exactarea)
return err
end

###############################################################
###############################################################
################################################################
################################################################

################################################################
################################################################
function greens3(nn,n)
#n points for each path
#n segments in x and y for ODE method to get integrand
#Actual integrand for this example is f(x,y)=x^2+y^2

n1=nn
n2=nn
n3=nn
n4=nn
#set parameters for ODE method
nx = n
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for ODE method
z = Array{Float64}(undef, nx, ny)

#euler's method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  z[j,i]=z[j,i-1]+hx*f(x1[i-1],y1[j])
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

#interpolate outer circle
t1 = 0:h1:1
itp1 = Interpolations.scale(interpolate(p1, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t1, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n1)
dy[1]=(p1[2,2]-p1[1,2])/h1
for i=2:n1
dy[i]=(p1[i,2]-p1[i-1,2])/h1
end
#interpolate dy
itpdy1 = linear_interpolation(t1, dy)
#approximate interal
int1, error = quadgk(t -> itp(itp1(t,1),itp1(t,2))*itpdy1(t), 0, 1)

#line integral for horizontal line
t2 = 0:h2:1
itp2 = Interpolations.scale(interpolate(p2, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t2, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n2)
dy[1]=(p2[2,2]-p2[1,2])/h2
for i=2:n2
dy[i]=(p2[i,2]-p2[i-1,2])/h2
end
#interpolate dy
itpdy2 = linear_interpolation(t2, dy)
#approximate interal
int2, error = quadgk(t -> itp(itp2(t,1),itp2(t,2))*itpdy2(t), 0, 1)

#line integral for inner circle
t3 = 0:h3:1
itp3 = Interpolations.scale(interpolate(p3, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t3, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n3)
dy[1]=(p3[2,2]-p3[1,2])/h3
for i=2:n3
dy[i]=(p3[i,2]-p3[i-1,2])/h3
end
#interpolate dy
itpdy3 = linear_interpolation(t3, dy)
#approximate interal
int3, error = quadgk(t -> itp(itp3(t,1),itp3(t,2))*itpdy3(t), 0, 1)



#line integral for line w/ slope 1
t4 = 0:h4:1
itp4 = Interpolations.scale(interpolate(p4, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t4, 1:2)
#approximate dy at knots
dy=Array{Float64}(undef,n4)
dy[1]=(p4[2,2]-p4[1,2])/h4
for i=2:n4
dy[i]=(p4[i,2]-p4[i-1,2])/h4
end
#interpolate dy
itpdy4 = linear_interpolation(t4, dy)
#approximate interal
int4, error = quadgk(t -> itp(itp4(t,1),itp4(t,2))*itpdy4(t), 0, 1)





int=(int1+int2+int3+int4)


exactarea=35*pi
err=abs((exactarea-int)/exactarea)
return err
end

###############################################################
###############################################################
################################################################
################################################################

function greens2(nn,n)
#nn points for each path
#n segments in x and y for ODE method to get integrand
#numeric approx of integral using Green's Thm
#Actual integrand for this example is f(x,y)=x^2+y^2

n1=nn
n2=nn
n3=nn
n4=nn
#set parameters for euler's method
nx = n
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for euler's method
z = Array{Float64}(undef, nx, ny)

#euler's method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  z[j,i]=z[j,i-1]+hx*f(x1[i-1],y1[j])
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

int1=0.0 #line integral for outer circle
for i=2:n1
x=p1[i,1]
y=p1[i,2]
dy=(p1[i,2]-p1[i-1,2])/h1
int1+=itp(x,y)*dy
end
int1*=h1

int2=0.0 #line integral for horizontal line
for i=2:n2
x=p2[i,1]
y=p2[i,2]
dy=(p2[i,2]-p2[i-1,2])/h2
int2+=itp(x,y)*dy
end
int2*=h2

int3=0.0 #line integral for inner circle
for i=2:n3
x=p3[i,1]
y=p3[i,2]
dy=(p3[i,2]-p3[i-1,2])/h3
int3+=itp(x,y)*dy
end
int3*=h3


int4=0.0 #line integral for line w/ slope 1
for i=2:n4
x=p4[i,1]
y=p4[i,2]
dy=(p4[i,2]-p4[i-1,2])/h4
int4+=itp(x,y)*dy
end
int4*=h4

int=(int1+int2+int3+int4)


exactarea=35*pi

err=abs((exactarea-int)/exactarea)
return err
end


################################################################
################################################################

function greens1(nn,n)  

#nn points for each path
#n segments in x and y to get integrand
#numeric approx of integral using Green's Thm 
#Actual integrand for this example is f(x,y)=x^2+y^2

n1=nn
n2=nn
n3=nn
n4=nn
#set parameters for ODE method
nx = n 
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for ODE method
z = Array{Float64}(undef, nx, ny)

#rk4 method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  k1=f(x1[i-1],y1[j])
  k2=f((x1[i-1]+x1[i])/2,y1[j])
  k4=f(x1[i],y1[j])
  z[j,i]=z[j,i-1]+(hx/6)*(k1+4*k2+k4)
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution
#itp = CubicSplineInterpolation((x1,y1), z; bc=Line(OnGrid()), extrapolation_bc=Throw())

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points 
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line 
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

#interpolate outer circle
t1 = 0:h1:1
###approximate dy at knots
dy=Array{Float64}(undef,n1)
dy[1]=(-3*p1[1,2]+4*p1[2,2]-p1[3,2])/(2*h1)
for i=2:n1-1
dy[i]=(p1[i+1,2]-p1[i-1,2])/(2*h1)
end
dy[n1]=(3*p1[n1,2]-4*p1[n1-1,2]+p1[n1-2,2])/(2*h1)
int1=itp(p1[1,1],p1[1,2]).*dy[1]+itp(p1[n1,1],p1[n1,2]).*dy[n1]
for i=2:n1-1
int1+=2*itp(p1[i,1],p1[i,2]).*dy[i]
end
int1*=(h1/2)

#line integral for horizontal line
t2 = 0:h2:1
#approximate dy at knots
dy=Array{Float64}(undef,n2)
dy[1]=(-3*p2[1,2]+4*p2[2,2]-p2[3,2])/(2*h2)
for i=2:n2-1
dy[i]=(p2[i+1,2]-p2[i-1,2])/(2*h2)
end
dy[n2]=(3*p2[n2,2]-4*p2[n2-1,2]+p2[n2-2,2])/(2*h2)
int2=itp(p2[1,1],p2[1,2]).*dy[1]+itp(p2[n2,1],p2[n2,2]).*dy[n2]
for i=2:n2-1
int2+=2*itp(p2[i,1],p2[i,2]).*dy[i]
end
int2*=(h2/2)

#line integral for inner circle
t3 = 0:h3:1
#approximate dy at knots
dy=Array{Float64}(undef,n3)
dy[1]=(-3*p3[1,2]+4*p3[2,2]-p3[3,2])/(2*h3)
for i=2:n3-1
dy[i]=(p3[i+1,2]-p3[i-1,2])/(2*h3)
end
dy[n3]=(3*p3[n3,2]-4*p3[n3-1,2]+p3[n3-2,2])/(2*h3)
int3=itp(p3[1,1],p3[1,2]).*dy[1]+itp(p3[n3,1],p3[n3,2]).*dy[n3]
for i=2:n3-1
int3+=2*itp(p3[i,1],p3[i,2]).*dy[i]
end
int3*=(h3/2)

#line integral for line w/ slope 1
t4 = 0:h4:1
#approximate dy at knots
dy=Array{Float64}(undef,n4)
dy[1]=(-3*p4[1,2]+4*p4[2,2]-p4[3,2])/(2*h4)
for i=2:n4-1
dy[i]=(p4[i+1,2]-p4[i-1,2])/(2*h4)
end
dy[n4]=(3*p4[n4,2]-4*p4[n4-1,2]+p4[n4-2,2])/(2*h4)
int4=itp(p4[1,1],p4[1,2]).*dy[1]+itp(p4[n4,1],p4[n4,2]).*dy[n4]
for i=2:n4-1
int4+=2*itp(p4[i,1],p4[i,2]).*dy[i]
end
int4*=(h4/2)

int=(int1+int2+int3+int4)


exactarea=35*pi

err=abs((exactarea-int)/exactarea)

return err

end

#####################################################
#####################################################

function greens_comp()
n=(10,20,40,80,160,320,640,1280)
n2=(5,10,20,40,80,160,320,640)
#####################################################
@printf "\n EXAMPLE 1: f(x,y)=x^2+y^2"
@printf "\n INTEGRAND:euler  FDM:1st order  QUADRATURE:left hand rule"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

er=Array{Float64}(undef,8)
eoc=Array{Float64}(undef,8)

er[1]=greens2(n[1],4*n[1])
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:8
er[i]=greens2(n[i],4*n[i])
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end

#####################################################
@printf "\n INTEGRAND:euler  FDM:1st order  QUADRATURE:quadgk"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

er=Array{Float64}(undef,8)
eoc=Array{Float64}(undef,8)

er[1]=greens3(n[1],4*n[1])
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:8
er[i]=greens3(n[i],4*n[i])
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end
#####################################################
@printf "\n INTEGRAND:rk4  FDM:1st order  QUADRATURE:quadgk"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

er=Array{Float64}(undef,8)
eoc=Array{Float64}(undef,8)

er[1]=greens4(n[1],4*n[1])
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:8
er[i]=greens4(n[i],4*n[i])
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end
#####################################################
@printf "\n INTEGRAND:rk4  FDM:2nd order  QUADRATURE:quadgk"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

er=Array{Float64}(undef,8)
eoc=Array{Float64}(undef,8)

er[1]=greens5(n[1],4*n[1])
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:8
er[i]=greens5(n[i],4*n[i])
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end


#####################################################
@printf "\n INTEGRAND:rk4  FDM:4th order  QUADRATURE:quadgk"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

er=Array{Float64}(undef,8)
eoc=Array{Float64}(undef,8)

er[1]=greens6(n[1],4*n[1])
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:8
er[i]=greens6(n[i],4*n[i])
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end



#####################################################
@printf "\n INTEGRAND:rk4  FDM:2nd order  QUADRATURE:trapezoidal rule"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

er=Array{Float64}(undef,8)
eoc=Array{Float64}(undef,8)

er[1]=greens1(n[1],4*n[1])
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:8
er[i]=greens1(n[i],4*n[i])
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end

#####################################################
@printf "\n INTEGRAND:rk4  FDM:4th order  QUADRATURE:trapezoidal rule"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

er=Array{Float64}(undef,8)
eoc=Array{Float64}(undef,8)

er[1]=greens7(n[1],4*n[1])
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:8
er[i]=greens7(n[i],4*n[i])
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end

end
