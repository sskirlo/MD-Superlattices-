

close all; 
clear all; 
clc; 

anb=3.3008; 
acu=3.615019;

Tnb=[(6^(1/2)-2)*(acu+anb)/anb, (3^(1/2)-2^(1/2))*(acu-2*anb)/2/anb; 
     (6-3*6^(1/2))/(3)^(1/2)+2*(3-6^(1/2))/3^(1/2)*acu/anb, (3-6^(1/2))+(1/2*6^(1/2)-1)*acu/anb        ];
 

Tcu=[-(6^(1/2)-3)*(acu+anb)/acu, (3^(1/2)-2^(1/2))*(acu-anb)/acu; 
     (6-2*6^(1/2))/(2)^(1/2)+3*(2-6^(1/2))/2^(1/2)*anb/acu,(6^(1/2)-2)-(6^(1/2)-3)*anb/acu              ];
 
Tcunb=[(3/2)^(1/2)*anb/acu, (2)^(1/2)/6*anb/acu ; 
          0,                     4/3*anb/acu      ]; 
      
%we want to check the action of the       
      
disp(Tnb*Tcunb); 
disp(Tcu); 
      

v1=acu/(2)^(1/2)*[1/2,3^(1/2)/2]; 
v2=acu/2^(1/2)*[1/2,-3^(1/2)/2];

v1p=[1/2,3^(1/2)/2]*(30-12*(6)^(1/2))^(1/2)*anb;
v2p=[(1/3)^(1/2),-(2/3)^(1/2)]*3/2*(6^(1/2)-2)*acu;

disp((Tcu*[v1.',v2.']).');
disp(v1p); 
disp(v2p); 

I=zeros(2,2); 
I(1,1)=1; 
I(2,2)=1; 

[V,D]=eig(I-Tcu); 
disp(V); 
disp(D); 


