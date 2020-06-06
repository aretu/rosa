function [coor1,period,amplitude,offset]=trapezoidal(slip,vmax,acc,fact)

% slip=1.5;
% vmax=3;
% acc=8;
ro=25/2000; ri=15/2000;
vmax_rpm=1/(4/3*pi*(ro^2+ro*ri+ri^2)/(ro+ri))*60*vmax;

base=slip/vmax-vmax/acc;
if or(base<0,vmax_rpm>1500)
    disp('error!!!!')
elseif base==0
    disp('triangular function!')
else
    Base=slip/vmax+vmax/acc;
    
    c=0.5*(Base-base);
    d=0.5*(Base-base)+base;
    
    A = [0 0] ;
    E = [Base*fact 0];
    B = [Base 0] ;
    C = [d vmax] ;
    D = [c vmax] ;
    coor = [A; E; B; C; D] ;
    figure(1000)
    patch(coor(:,1), coor(:,2),'-or'); xlabel('Time (s)'); ylabel('Velocity (m/s)');
    
    A1=[0 0];
    E1=[1 0];
    B1=[Base/Base/fact 0];
    C1=[d/Base/fact vmax/vmax];
    D1=[c/Base/fact vmax/vmax];
    coor1 = [A1 ; E1; B1; C1; D1] ;
    figure(999)
    patch(coor1(:,1), coor1(:,2),'-og'); xlabel('Norm units'); ylabel('Norm units');
    
    
    period=Base; period=period*fact;
    amplitude=vmax_rpm/149.9;
    offset=amplitude/2;
end
end