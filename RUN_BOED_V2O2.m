clear
clc

% 1) DEFINE INPUT SPACE 
VR1 = input('input the array of possible values of vaiable 1: \n');  %Variable 1
VR2 = input('input the array of possible values of vaiable 2: \n'); % Variable 2
% 2) DEFINE OBJECTIVE FUNCTIONS
objectfun = input('input the function handles of two objectives: \n'); % Define Objective functions
% 3) DEFINE nE and nI
nE= input('input the number of evalutions of the BOED process: \n'); % Preselected number of evaluations of he BOED process. Each evaluation correspods to an experiment
nI= input('input the initial number: \n'); % Number of randomly selected experiments
reference_point = input('input the reference point for the objectives domain: \n');


% CODE BEGINNING 
XSPACE=DXSPACE(VR1,VR2); % DEFINE XSPACE MATRIX BASED ON THE INPUT SPACE
[X, Y, XSPACE, qnum, dnum] = DEVELOP_INITIAL_IO_DATABASE(objectfun,XSPACE,nI); % DEVELOP INITIAL INPUT-OUTPUT DATABASE
[X, Y, XSPACE, PARETO_FRONT, qnum, dnum]= MAIN_BOED_ITERATIVE_LOOP(objectfun,XSPACE,X,Y,reference_point,nI,nE,qnum, dnum); % MAIN BOED ITERATIVE LOOP
PLOT_FIGURES(X,Y,XSPACE,PARETO_FRONT,qnum, dnum)


%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################      
%######################################################################################  MAIN BOED FUNCTIONS    ###############################################################################################
%##############################################################################################################################################################################################################      
%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################  


function [X, Y, XSPACE, qnum, dnum] = DEVELOP_INITIAL_IO_DATABASE(objectfun,XSPACE,nI)
%qnum the number of quried rows of YSPACE 
%dnum the number of discarded data 888 999
 sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  DEVELOP INITIAL INPUT-OUTPUT DATABASE  ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])    

%INITIALIZE VARIABLES AND MATRICES
qnum = 0;
dnum = 0;
Y = nan(1, 2); %INITIALIZE YSPACE MATRIX
X = nan(1, 2); %INITIALIZE XSPACE MATRIX

m=1;     
while m<=nI 
            
    UNEXPLORED_XSPACE=XSPACE(qnum+1:size(XSPACE,1)-dnum,2:3);
    n=size(UNEXPLORED_XSPACE,1); % NUMBER OF REMAINING POTENTIAL EXPERIMENTS
    
    if (n==0) % if m reaches the size of XSPACE -  the number of discarded data +1 break the while loop
    sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  NO MORE INITIAL EXPERIMENTS CAN BE PERFORMED ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])    
    break
    end   
    clear data
    
    % SELECTION OF INPUT VARIABLES WITHOUT REPLACEMENT FROM THE PREDIFINED INPUT SPACE (XSPACE)
    INPUT=datasample(UNEXPLORED_XSPACE,1);
       
    % CALCULATION OF OPERATIONAL OBJECTIVES
    NI_COMPO=INPUT(1,1);
    VF=INPUT(1,2);
    [idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(NI_COMPO,VF,XSPACE);
    [OBJECTIVES_VALUE, ERROR]=CALCULATE_OBJECTIVES(objectfun,NI_COMPO,VF);
    
    % UPDATE INPUT-OUTPUT DATABASE
    [X,Y]=UPDATE_IO_DATABASE(OBJECTIVES_VALUE,X,Y,INPUT,ERROR);
    [XSPACE, dnum, qnum ] =SORT_XSPACE(XSPACE,dnum,qnum,idx,ERROR);                                            
    
    if ERROR==0
    sprintf(['RANDOMLY SELECTED EXPERIMENT: ',num2str(m),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' SUCCESFULLY COMPLETED'])
    m=m+1;
    else
    sprintf(['RANDOMLY SELECTED EXPERIMENT: ',num2str(m),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' NOT COMPLETED - PROBLEMATIC INPUT VALUE. SELECT NEXT MATERIAL TO TEST FROM THE REDUCED INPUT SPACE'])    
    end
        
end   

sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  INITIAL_DATABASE COMPLETED ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])

end



function [X, Y] =UPDATE_IO_DATABASE(OBJECTIVES_VALUE,X,Y,INPUT,ERROR)

m=size(X,1);

if ERROR==0 % UPDATE ONLY WHEN THERE IS NO ERROR
    
    if m==1 && isnan(X(1,1))==1 % For the first time that X,Y are updated
        X(1,1:2)=INPUT(1:2);
        Y(1, 1:2) = OBJECTIVES_VALUE;
    elseif m>=1 && isnan(X(1,1))==0
        X(m+1,1:2)=INPUT(1:2);
        Y(m+1, 1:2) = OBJECTIVES_VALUE;
    end

end
   

end

function [X, Y,XSPACE, PARETO_FRONT, qnum, dnum]=MAIN_BOED_ITERATIVE_LOOP(objectfun,XSPACE,X,Y,boundpoint,nI,nE,qnum, dnum)
sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  BOED MAIN ITERATIVE LOOP  ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])

nT=size(XSPACE,1);
nRT=nT-nI-dnum; %  Remaining INPUT SPACE that can be covered during the evaluations of the BOED method
i=1;
 
while i<=nE  && i <=nRT    % IT WILL STOP WHEN THE STOPNUM IS REACHED OR WHEN THE ENTIRE SPACE IS COVERED. THE ENTIRE SPACE IS DEFINED AS THE INITIAL_SPACE-THE INITIAL RUNS-NOT ACCEPTABLE POINTS
   
         
    [Y_MEAN, Y_STD]= MACHINE_LEARNING(X,Y,XSPACE, qnum, dnum); % Construct Surrogate model based on the input-output data (X,Y)
    [INPUT , ~]= SELECTOR(X,Y, Y_MEAN, Y_STD,XSPACE, qnum, dnum, boundpoint); % Select the next material to test (INPUT)
    NI_COMPO=INPUT(1,1);
    VF=INPUT(1,2);
    [idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(NI_COMPO,VF,XSPACE);
    [OBJECTIVES_VALUE, ERROR]=CALCULATE_OBJECTIVES(objectfun,NI_COMPO,VF);
    
    % UPDATE INPUT-OUTPUT DATABASE
    [X,Y]=UPDATE_IO_DATABASE(OBJECTIVES_VALUE,X,Y,INPUT,ERROR);
    [XSPACE, dnum, qnum ] =SORT_XSPACE(XSPACE,dnum,qnum,idx,ERROR);   
    
    % UPDATE AUXILARY VARIABLES
     nRT=nT-nI-dnum;
    
    if ERROR==0
    sprintf(['EVALUATION: ',num2str(i),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' SUCCESFULLY COMPLETED'])
    i=i+1;
    else
    sprintf(['EVALUATION ',num2str(i),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' PROBLEMATIC INPUT VALUE, PROCEED TO NEXT CASE_NUMBER'])    
    end
    
end
sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  BOED PROCESS COMPLETED  ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])
PARETO_FRONT = FIND_PARETO_FRONT(Y);

end

function [Y_MEAN, Y_STD]= MACHINE_LEARNING(X,Y,XSPACE, qnum, dnum)
XSPACER=XSPACE(:,2:3);

MEAN=mean(XSPACER);
STD=std(XSPACER);

X_NORM = (X-MEAN)./(ones(size(X, 1),1)*STD);
left_idx = qnum+1:size(XSPACER, 1)-dnum;% from current number to the end of xspace
X_TEST=XSPACER(left_idx, :);
X_TEST_NORM=(X_TEST-MEAN)./(ones(size(X_TEST, 1),1)*STD);
[Y_MEAN, Y_STD] = GPR_SURROGATE_MODEL(X_NORM,Y,X_TEST_NORM);
end

function [Y_MEAN, Y_STD] = GPR_SURROGATE_MODEL(X, Y, X_TEST)
% ymean,  size R*C 
% ysd,  Standard deviation, sigma. ysd^2 = variance 
C = size(Y, 2);
R = size(X_TEST, 1);
Y_MEAN = zeros(R, C);
Y_STD = zeros(R, C);
for n = 1:C
    gprMd = fitrgp(X, Y(:,n));%, 'BasisFunction', 'linear');%
    for m = 1:R
        [Y_MEAN(m, n), Y_STD(m, n), ~] = predict(gprMd, X_TEST(m, :));% 
    end
end
end

function [INPUT, p_integral]= SELECTOR(~,Y,Y_MEAN, Y_STD,XSPACE, qnum, dnum, boundpoint)
XSPACER=XSPACE(:,2:3);
left_idx = qnum+1:size(XSPACER, 1)-dnum;% from current number to the end of xspace
X_TEST=XSPACER(left_idx, :);
m=size(left_idx,2);
PARETO_FRONT = FIND_PARETO_FRONT(Y);% value of y object, sorted
max_integral = 0;

for i = 1:m
    
            
    p_integral = EHVI(Y_MEAN(i,:), Y_STD(i,:), PARETO_FRONT, boundpoint);
    
    if p_integral >= max_integral
        max_integral = p_integral;
        NI_COMPO=X_TEST(i,1);
        VF=X_TEST(i,2);
        [idx, ~]=IDENTIFY_IDX_AND_CASE_NUMBER(NI_COMPO,VF,XSPACE);
        next_idx = idx;
    end
    
end

INPUT=XSPACER(next_idx,:);


end

function p_integral = EHVI(mu, sigma, PARETO_FRONT, extreme_point)

Phi = @(s) 1/2*(1+erf(s/sqrt(2)));
bigpsi = @(a, b, mu, sigma) sigma*normpdf((b-mu)/sigma)+(a-mu)*Phi((b-mu)/sigma);
n = size(PARETO_FRONT, 1);

part1 = zeros(1, n);
part2 = zeros(1, n);

y = zeros(n+2, 2);
y(1, :) = [extreme_point(1), -inf];
y(2:end-1, :) = PARETO_FRONT;
y(end, :) = [-inf, extreme_point(2)];

% y is the objective front array
for i = 2:size(y, 1)
    part1(i-1) = Phi((y(i,1)-mu(1))/sigma(1))*(y(i-1, 1)-y(i, 1))*...
        bigpsi(y(i, 2), y(i,2), mu(2), sigma(2));
    part2(i-1) = ( bigpsi(y(i-1, 1), y(i-1, 1), mu(1), sigma(1)) - ...
        bigpsi(y(i-1, 1), y(i, 1), mu(1), sigma(1)) )*...
        bigpsi(y(i, 2), y(i, 2), mu(2), sigma(2));
end
part1(end) = 0;
p_integral = sum(part1)+sum(part2);

end


%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################      
%#################################################################################     AUXILARRY FUNCTIONS         ############################################################################################
%##############################################################################################################################################################################################################      
%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################   

function [idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(VR1,VR2,XSPACE)
% IDENTIFIES WHICH ROW OF XSPACE REPRESENTS THE INPUT VARIABLES VR1,VR2
% IT RETURNS THE CASE_NUMBER THAT CORRESPONDS TO THE IDENTIFIED ROW
% DONT CONFUSE IDX WITH THE CASE_NUMBER.

idx = 0;
for m = 1:size(XSPACE, 1)
    if XSPACE(m, 2) == VR1 && XSPACE(m, 3) == VR2 
        idx = m;
        CASE_NUMBER=XSPACE(m, 1);
    end
end
  

if idx == 0
error('something wrong')
end


end

function XSPACE=DXSPACE(A1,A2)
delete('CASES.xlsx');
xspace1=COMB(A1,A2);
XSPACE(:,1)=1:1:size(xspace1,1); % THIS COLUMN NOTES THE CASE NUMBERS.
XSPACE(:,2)=xspace1(:,2);
XSPACE(:,3)=xspace1(:,1);
XSPACE(:,2)=round(XSPACE(:,2),2);
XSPACE(:,3)=round(XSPACE(:,3),3);
xlswrite('CASES.xlsx',XSPACE);

end

function A=COMB(A1,A2)
% TAKES AS INPUT TWO MATRIXES A1 AND A2 AND RETURNS THE COMBINATIONS OF
% THEM

if size(A1,1)==1 % WHEN THE A1 IT IS VECTOR AND IT IS DEFIND AS A ROW I HAVE TO SWITCH IT TO COLUMN
A1=A1';
end


index1=max(size(A1));
index2=min(size(A1));
index3=max(size(A2));
index4=min(size(A2));


k=0;
for i =1:index1
    for j=1:index3
    start=(j-1)*index1+1;
    finish=(j)*index1;    
    A(start:finish,1:index2)=A1(:,1:index2);
    A(start:finish,index2+1)=A2(j);
    k=k+1;
    end
end


end

function [ XSPACE, dnum, qnum ] =SORT_XSPACE(XSPACE,dnum,qnum,idx,ERROR)

    if ERROR~=0
        XSPACE([end-dnum, idx], :) = XSPACE([idx, end-dnum], :);
        dnum = dnum+1;
    else
        qnum = qnum+1;
        XSPACE([qnum, idx], :) = XSPACE([idx, qnum], :);
    end
   
end

function output_array = FIND_PARETO_FRONT(input_array)
temp = 1:size(input_array, 1);
for y_idx = 1:size(input_array, 1)
    for y_idx2 = 1:size(input_array, 1)
        if input_array(y_idx, :)>input_array(y_idx2, :)
            % if it dominate other point, it's not small enough, should be
            % excluded
            temp = setdiff(temp, y_idx);
        end
    end
end
output_array = input_array(temp, :);
[~, input_idx] = sort(output_array(:,2));
output_array = output_array(input_idx, :);
end

function PLOT_FIGURES(X,Y,~,PARETO_FRONT,qnum, ~)

OBJECTIVE1_max=max(Y(:,1));
OBJECTIVE1_min=min(Y(:,1));
OBJECTIVE2_max=max(Y(:,2));
OBJECTIVE2_min=min(Y(:,2));
scalex=1.24;
scaley=1.15;


figure(1)
% PLOT FORMATING
x0=0;
y0=0;
width=500;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
set(gca,'fontsize',28)
set(gca,'fontWeight','bold')
set(gca,'fontName','Arial')
set(gca,'linewidth',2)
set(gca, 'box', 'on')
points=500;% THE NUMBER OF THE POINTS THAT I WILL DISCRITIZE THE CMAP
xlabel('Objective 1','FontSize', 36)
ylabel('Objective 2','FontSize', 36)
axis([0 OBJECTIVE1_max*scalex 0 OBJECTIVE2_max*scaley])
daspect([OBJECTIVE1_max*scalex OBJECTIVE2_max*scaley 1])

% DEFINITIO OF THE VARIABLES
hold on
NI_COMPO_max=max(X(:,1));
NI_COMPO_min=min(X(:,1));
VF_max=max(X(:,2));
VF_min=min(X(:,2));


% STORE EVERYTHING IN THE OUT VARIABLES FOR DEBUGGING
clear out
for i=1:qnum
NI_COMPO=X(i,1);
VF=X(i,2);
color=lin(NI_COMPO,NI_COMPO_max,NI_COMPO_min,1,0);
rad=lin(VF,VF_max,VF_min,80,20);    
out(i,1)=Y(i,1);
out(i,2)=Y(i,2);
out(i,3)=VF;
out(i,4)=NI_COMPO;
out(i,5)=fin(color,points);
end

% DEFINITION OF THE COLOR MAP
clear cmap
cmap = jet(points);
%cmap=flip(cmap)
c=colormap(cmap);
caxis([NI_COMPO_min NI_COMPO_max]);
c=colorbar;
c.Label.String = 'Ni Concentration (% at.)';
c.Label.FontSize = 36;
c.Label.VerticalAlignment='top';

% PLOT ENTIRE OBJECTIVE SPACE
for i=1:qnum
NI_COMPO=X(i,1);
VF=X(i,2);
color=lin(NI_COMPO,NI_COMPO_max,NI_COMPO_min,1,0);
rad=lin(VF,VF_max,VF_min,90,30);    
hold on
plot(Y(i,1), Y(i,2),'MarkerSize', rad,'marker','.','Color',cmap(fin(color,points),:))% Show the objectives values for all the executed analysis
hold on
plot(Y(i,1), Y(i,2),'MarkerSize', rad/3.3,'marker','o','Color','black')% Show the objectives values for all the executed analysis
end


%PLOT PARETO FRONT
for i=1:size(PARETO_FRONT,1)
    
for j=1:qnum
    if PARETO_FRONT(i,1)==Y(j,1) && PARETO_FRONT(i,2)==Y(j,2)
    index=j;
    end
end
NI_COMPO=X(index,1);
VF=X(index,2);
color=lin(NI_COMPO,NI_COMPO_max,NI_COMPO_min,1,0);
rad=lin(VF,VF_max,VF_min,90,30);
hold on
plot(PARETO_FRONT(i,1), PARETO_FRONT(i,2)); %  Show Front Array
end
plot(PARETO_FRONT(:,1), PARETO_FRONT(:,2),'-.black','LineWidth', 2); %  Show Front Array

end

function value=lin(X_EVAL,X_MAX,X_MIN,Y_MAX,Y_MIN)
% DEFINE A LINEAR FUNCTION WHICH AT X_MAX -->Y_MAX AND X_MIN-->Y_MIN
% CALCULATE THE VALUE OF THE FUNCTION AT X_EVAL
   alpha=(Y_MAX-Y_MIN)/(X_MAX-X_MIN);
   bhta=Y_MIN-alpha*X_MIN;
   value=alpha*X_EVAL+bhta;

end

function output=fin(INPUT,L)
%DEFINES AN XX SPACE FROM 0 TO 1 WITH POINTS EQUAL TO THE ONES SPECIFIED BY L.
%THEN IT CHECKS THE DEFINED INPUT BETWEEN WHICH POINTS BELONGS IN THE XX
%SPACE AND IF FOR EXAMPLE IT BELONGS BETWEEN 4 AND 5 IT RETURNS THE 4 (THE FIRST OF THE TWO).
% THIS IS USEFULL FOR THE DEFINITION OF THE COLORS. THE CMAP IS DEFINED AS
% A COLOR MAP WITH L POINTS. THE 1ST POINT IS THE LOW THE LAST IS THE HIGH.
% HENCE FIRST I MAP THE COMPOSITION FROM 1 TO 0 AND THEN I DEFINE USING THE
% FIN THE COLOR THAT CORRESPONDS TO THAT COMPOSITION. FOR EXAMPLE IF
% COMPOSITION 50.3 CORRESPONDS TO COLOR 0.64. THEN I SEARCH IN THE XX TO
% FIND BETWEEN WHICH POINTS THE 0.64 LIES AND THEN I PICK THE INDEX OF THE
% FIRST OF THE TWO. THEN I DEFINE AS COLOR THE COLOR FROM THE CMAP THAT
% CORRESPONDS TO THAT INDEX.


xx=0:1/(L-1):1;

for i=1:L

if (INPUT>=xx(i))  
  output=i;  
end
    
end

end


