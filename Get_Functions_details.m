

function [lu,fobj,threshold] = Get_Functions_details(F, dim)


switch F
    case 'F1'
        fobj = @F1;
        
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 0.1;
        
        
    case 'F2'
        fobj = @F2;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 0.1;
        
        
    case 'F3'
        fobj = @F3;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 0.1;
        
    case 'F4'
        fobj = @F4;
        lu = [-10 * ones(1, dim); 10 * ones(1, dim)];
        threshold = 0.1;
        
    case 'F5'
        fobj = @F5;
        lu = [-30 * ones(1, dim); 30 * ones(1, dim)];
        threshold = 100;
        
    case 'F6'
        fobj = @F6;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 0.1;
    case 'F7'
        fobj = @F7;
        lu = [-10 * ones(1, dim); 10 * ones(1, dim)];
        threshold = 100;
    case 'F8'
        fobj = @F8;
        lu = [-5.12 * ones(1, dim); 5.12 * ones(1, dim)];
        threshold = 50;
    case 'F9'
        fobj = @F9;
        lu = [-5 * ones(1, dim); 10 * ones(1, dim)];
        threshold = 0.1;
    case 'F10'
        fobj = @F10;
        lu = [-10 * ones(1, dim); 10 * ones(1, dim)];
        threshold = 0.1;
    case 'F11'
        fobj = @F11;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 100;
    case 'F12'
        fobj = @F12;
        lu = [-32 * ones(1, dim); 32 * ones(1, dim)];
        threshold = 0.1;
    case 'F13'
        fobj = @F13;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 1;
    case 'F14'
        fobj = @F14;
        lu = [-600 * ones(1, dim); 600 * ones(1, dim)];
        threshold = 0.1;
    case 'F15'
        fobj = @F15;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 0.1;
    case 'F16'
        fobj = @F16;
        lu = [-10 * ones(1, dim); 10 * ones(1, dim)];
        threshold = 100;
    case 'F17'
        fobj = @F17;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 1;
    case 'F18'
        fobj = @F18;
        lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
        threshold = 1;
        
end
lu = lu';
end

% F1 Sphere Function

function o = F1(x)
o=sum(x.^2,2);
end

% F2 Schwefel 1.2 Function

function o = F2(x)
dim=size(x,2);
swarm=size(x,1);
o=zeros(swarm,1);
for i=1:dim
    o=o+sum(x(:,1:i),2).^2;
end
end


% F3 Schwefel 2.21 Function

function o = F3(x)
dim=size(x,2);
swarm=size(x,1);
o=zeros(swarm,1);
for i = 1:swarm
    o(i)=max(abs(x(i,1:dim)));
end
end

% F4 Schwefel 2.22 Function

function o = F4(x)
o=sum(abs(x),2)+prod(abs(x),2);
end

% F5 Rosenbrock’s Function

function o = F5(x)
dim=size(x,2);
o=100*sum((x(:,2:dim)-(x(:,1:dim-1)).^2).^2,2)+sum((x(:,1:(dim-1))-1).^2,2);

end

% F6 Step Function

function o = F6(x)
dim=size(x,2);
o=sum(floor((x(:,1:dim)+.5)).^2,2);
end

% F7 Bent Cigar Function
function o = F7(x)
dim=size(x,2);
o=x(1)^2+10^6*(sum(x(:,2:dim).^2,2));
end

% F8 Rastrigin’s Function
function o = F8(x)
dim=size(x,2);
o=sum(x(:, 1:dim).^ 2 - 10 .* cos(2 .* pi .*x(:,1:dim)) + 10, 2);
end

% F9 Zakharov Function
function o = F9(x)
dim=size(x,2);
sum_1=0;
for i=1:dim
    sum_1=sum(0.5*i*x(i),2);
end
o=sum(x.^2,2)+sum_1^2+sum_1^4;
end

% F10 Levy Function
function o = F10(x)
dim=size(x,2);
w=1+(x-1)/4;
sum_1=sum(((w(:,1:dim-1)-1).^2).*(1+(10*((sin(pi*w(:,1:dim-1)+1)).^2))),2);
sum_2=(sin(pi*w(1)))^2;
sum_3=((w(dim)-1)^2)*(1+10*((sin(pi*w(dim)))^2));
o=sum_1+sum_2+sum_3;
end

% F11 High Conditioned Elliptic Function
function o = F11(x)
dim=size(x,2);
o=sum(((10^6).^(((1:dim)-1)/(dim-1))).*x(:,1:dim).^2,2);
end

% F12 Ackley’s Function
function o = F12(x)
dim=size(x,2);  
o=20+exp(1)-20*exp(-0.2*(1/dim*sum(x.^2,2)).^0.5)-exp(1/dim*sum(cos(2*pi*x),2));
end

% F13 HappyCat Function
function o = F13(x)
dim=size(x,2);
sum_1=0;
sum_2=0;
temp1=0;
temp2=0;
for i=1:dim
sum_1=sum_1+x(i)^2;
sum_2= sum_2+x(i);
temp1 = (abs(sum_1-dim)).^(1/4);
temp2 = (0.5* sum_1+ sum_2)/dim;
o(i)= temp1+ temp2+0.5;
end
end

% F14 Griewank’s Function
function o = F14(x)
dim=size(x,2);
swarm=size(x,1);
for i=1:dim
o=sum(x.^2,2)/4000-prod(cos(x./(repmat(1:dim,swarm,1)).^0.5),2)+1;
end
end
% F15 Sum of Different Power Function
function o = F15(x)
dim=size(x,2);
for i=1:dim
  o=sum(abs(x(:,i)).^(i+1),2);
end
end

% F16 Discus Function
function o = F16(x)
dim=size(x,2);
tem1=10.^6*((x(1)).^2);
tem2=sum ((x(:,2:dim)).^2,2);
o=tem1+tem2;
end

% F17 HGBat Function
function o = F17(x)
dim=size(x,2);
sum_1 = 0;
sum_2=0;
for i=1:dim
sum_1=sum_1+x(i)^2;
sum_2= sum_2+x(i);
o(i)=(abs(sum_1.^2- sum_2.^2)).^(1/2)+(0.5* sum_1+ sum_2)/dim+0.5;
end
end

% F18 Schaffer's F7 Function
function o = F18(x)
dim=size(x,2);
o=(sum((x(:,1:dim-1).^2+x(:,2:dim).^2).^0.25+((x(:,1:dim-1).^2+x(:,2:dim).^2).^0.25).*(sin(50*((x(:,1:dim-1).^2+x(:,2:dim).^2).^0.1))).^2,2))/(dim-1);    
end
