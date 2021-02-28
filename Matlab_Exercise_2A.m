a=input('Choose a question: ','s');
switch a
    %robot library of Corke MATLAB Robotics Toolbox must be installed to
    %use some functions 
case 'a'
    %Question A
    %Creating a function that returns a rotation matrix
    
    i1=[10  20 30]; %[alpha,beta,gamma]
    i2=[30 90 -55];
    
    RB_A_i1=rotate(i1); %rotation matrix i
    RB_A_i2=rotate(i2); %rotation matrix j
    
    %To demonstrate the six contraints for unitary orthonormal rotation
    %matricies , we use equation 2.59
    
    
    X1=RB_A_i1(:,1);
    Y1=RB_A_i1(:,2);
    Z1=RB_A_i1(:,3);
    
    check1=(norm(X1)==1);
    check2=(norm(Y1)==1);
    check3=(norm(Z1)==1);
    check4=(dot(X1,Y1)<10^-10);
    check5=(dot(X1,Z1)<10^-10);
    check6=(dot(Y1,Z1)<10^-10);
    
    if check1 & check2 & check3 & check4 & check5 & check6
        disp("All constraints are achieved");
    else
        disp("Not all constraints are achieved");
    end
    
    %To demonstrate the beautiful property , we check if RAinB =
    %RBinAinv=transpose(RBinA)
    A=transpose(RB_A_i1);
    B=inv(RB_A_i1);
    
    if (round(A,4)==round(B,4))
        disp("transpose(RA_B)==inv(RA_B) ==> RB_A == transpose(RA_B)==inv(RA_B)");
    end
case 'b'
    %Question B
    
    %here we will use the function antirotate on resulting matrix of i and ii
    
    [R1,R2]=(antirotate(RB_A_i1));
    disp("The first solution of i is ");disp(R1);
    disp("The second solution of i is ");disp(R2);
    
    %using R1 and R2 again in A we should get RB_A_i1
    
    if rotate(R1)-RB_A_i1<10^-10
        disp("First solution of i is verified")
    else
        disp("First solution of i is not verified")
    end
    
    if rotate(R2)-RB_A_i1<10^-10
        disp("Second solution of i is also correct")
    else
        disp("Second solution of i is also not correct")
    end
    
    % [R1,R2]=(antirotate(RB_A_i2));
    % disp("The first solution of ii is ");disp(R1);
    % disp("The second solution of ii is ");disp(R2);
    % in case beta=+-90 , we could only find the sum or difference of the
    % other two angles
case 'c'
    %Question 3
    
    i3=[0 20 0];
    R3=rotate(i3);
    P_B=[1 0 1]';
    A_P=R3*P_B %using rule 2.96
    plotter(P_B);
    plotter(A_P);
    grid on;
 case  'd'
     %Question 4
     
     V_R1=rpy2tr(i1,'deg','zyx'); % verify R1
     V_R2=rpy2tr(i2,'deg','zyx'); % verify R2
     
     if norm(RB_A_i1-V_R1(1:3, 1:3))<10^-10
         disp("The rotate algorithm is verified for i");
     else
         disp("The rotate algorithm is not verified for i");
     end
     
     if norm(RB_A_i2-V_R2(1:3, 1:3))<10^-10
         disp("The rotate algorithm is verified for ii");
     else
         disp("The rotate algorithm is not verified for ii");
     end
     
     V_i1=tr2rpy(RB_A_i1,'deg','zyx');
     V_i2=tr2rpy(RB_A_i2,'deg','zyx');
     
     if norm(V_i1-i1)<10^-10
         disp("The antirotate algorithm is verified for i");
     else
         disp("The antirotate algorithm is not verified for i");
     end
     
     if norm(V_i2-i2)<10^-10
         disp("The antirotate algorithm is verified for ii");
     else
         disp("The antirotate algorithm is not verified for ii");
     end
     
     %to verify c ,we use roty(beta)
     
     A_P_V=roty(20,'deg')*P_B;
     
     if norm(A_P_V-A_P)<10^-10
         disp("The rotate algorithm is verified for c");
     else
         disp("The rotate algorithm is not verified for c");
     end

end
function [R] = rotate(i)
   a=i(1,1);
   b=i(1,2);
   c=i(1,3);
   
   r11=cosd(a)*cosd(b);
   r12=cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c);
   r13=cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c);
   r21=sind(a)*cosd(b);
   r22=sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c);
   r23=sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c);
   r31=-sind(b);
   r32=cosd(b)*sind(c);
   r33=cosd(b)*cosd(c);
   
   R=[r11 r12 r13;r21 r22 r23;r31 r32 r33];
end

function [R1 , R2] = antirotate(RB_A)
% here we use the equation 2.66 to calculate the 2 solutions 
% we have a risk of having cos(b) =0 thus we solve this by isolating this
% case , thus having beta=+-90 will cause first column to be [0 0 +-1] and
% last row to be [+-1;0;0] , it is enough to compare the first column
if isequal(RB_A(:,1),[0;0;1]) 
    R1=[0 -90 -atan2d(RB_A(1,2),RB_A(2,2))];
    R2=R1;
    return;
end

if isequal(RB_A(:,1),[0;0;-1]) 
    R1=[0 90 atan2d(RB_A(1,2),RB_A(2,3))];
    R2=R1;
    return;
end
% the other solution would use bnew=pi-b

b=atan2d(-RB_A(3,1),sqrt(RB_A(1,1)^2+RB_A(2,1)^2));
a=atan2d(RB_A(2,1)/cosd(b),RB_A(1,1)/cosd(b));
c=atan2d(RB_A(3,2)/cosd(b),RB_A(3,3)/cosd(b));
R1=[a b c];

bnew=180-b;
anew=atan2d(RB_A(2,1)/cosd(bnew),RB_A(1,1)/cosd(bnew));
cnew=atan2d(RB_A(3,2)/cosd(bnew),RB_A(3,3)/cosd(bnew));
R2=[anew bnew cnew];
end

function [] = plotter(R)

plot3([0 R(1)],[0 R(2)],[0 R(3)]);
hold on;
end
