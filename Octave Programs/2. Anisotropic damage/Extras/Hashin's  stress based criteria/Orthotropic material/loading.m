function load=loading(ltype,dt,t,lam)

if ltype==1
    xa=t(1):dt:t(2);
    xa1=t(1); ya1=lam(1);
    xa2=t(2); ya2=lam(2);
    ya=(ya2-ya1)/(xa2-xa1)*(xa-xa1)+ya1;
    load=ya;    
elseif ltype==2
    xa=t(1):dt:t(2);
    xa1=t(1); ya1=lam(1);
    xa2=t(2); ya2=lam(2);
    ya=(ya2-ya1)/(xa2-xa1)*(xa-xa1)+ya1;

    xb=(t(2)+dt):dt:t(3);
    xb1=t(2); yb1=lam(2);
    xb2=t(3); yb2=lam(1);
    yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1;

    load=[ya yb];
    
elseif ltype==3
    xa=t(1):dt:t(2);
    xa1=t(1); ya1=lam(1);
    xa2=t(2); ya2=lam(2);
    ya=(ya2-ya1)/(xa2-xa1)*(xa-xa1)+ya1;

    xb=(t(2)+dt):dt:t(3);
    xb1=t(2); yb1=lam(2);
    xb2=t(3); yb2=lam(3);
    yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1;
    
    load=[ya yb];

elseif ltype==4
    xa=t(1):dt:t(2);
    xa1=t(1); ya1=lam(1);
    xa2=t(2); ya2=lam(2);
    ya=(ya2-ya1)/(xa2-xa1)*(xa-xa1)+ya1;

    xb=(t(2)+dt):dt:t(3);
    xb1=t(2); yb1=lam(2);
    xb2=t(3); yb2=lam(3);
    yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1;

    xc=(t(3)+dt):dt:t(4);
    xc1=t(3); yc1=lam(3);
    xc2=t(4); yc2=lam(1);
    yc=(yc2-yc1)/(xc2-xc1)*(xc-xc1)+yc1;

    load=[ya yb yc];
    
elseif ltype==5
    xa=t(1):dt:t(2);
    xa1=t(1); ya1=lam(1);
    xa2=t(2); ya2=lam(2);
    ya=(ya2-ya1)/(xa2-xa1)*(xa-xa1)+ya1;

    xb=(t(2)+dt):dt:t(3);
    xb1=t(2); yb1=lam(2);
    xb2=t(3); yb2=lam(3);
    yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1;

    xc=(t(3)+dt):dt:t(4);
    xc1=t(3); yc1=lam(3);
    xc2=t(4); yc2=lam(2);
    yc=(yc2-yc1)/(xc2-xc1)*(xc-xc1)+yc1;
    
    xd=(t(4)+dt):dt:t(5);
    xd1=t(4); yd1=lam(2);
    xd2=t(5); yd2=lam(3);
    yd=(yd2-yd1)/(xd2-xd1)*(xd-xd1)+yd1;
    
    xe=(t(5)+dt):dt:t(6);
    xe1=t(5); ye1=lam(3);
    xe2=t(6); ye2=lam(2);
    ye=(ye2-ye1)/(xe2-xe1)*(xe-xe1)+ye1;
    
    xf=(t(6)+dt):dt:t(7);
    xf1=t(6); yf1=lam(2);
    xf2=t(7); yf2=lam(1);
    yf=(yf2-yf1)/(xf2-xf1)*(xf-xf1)+yf1;

    load=[ya yb yc yd ye yf];
    
end
