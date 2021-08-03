% EPICS=Effective Pairwise Interactions for predicting Community Structure

%The following function takes n(number of species in the community),
%species abundances in monoculture and leave-one-out subcommunities and
%produces effective pairwise interactions and community structure using EPICS as
%outputs. In this code, we assumed that the microbial community follows
%the generalized Lotka-Volterra model. Other population dynamics models can be
%used instead.

clear all;
close all;
clc;
disp('Welcome to EPICS...');

%number of species in the community
number_of_species=input('\nInput number of species in the community=');


%monoculture abundances
ind_abun=input('\nInput monoculture abundances \n(Put your entries wrapped in a square bracket. \nPut a space or a comma after each entry)=');

%leave-one-out abundances
loo_abun=input('\nInput abundances in leave-one-out subcommunities \n(Rows are subcommunities, columns are species)=');

global n %number of species in the community


n=number_of_species;
for i=1:n
    A(i,i)=-1/ind_abun(i);%the diagonal elements are determined from monoculture abundances
    diagonalA(i)=A(i,i);
end

r=ones(n,1);
%This function computes effective pairwise interactions using EPICS

%intrinsic growth rate vector
r=1*ones(n,1)/1;

for i=1:n
self_intn(i,1)=-1/ind_abun(1,i);%the diagonal elements or self interaction terms are determined from monoculture abundances
end

% To obtain effective pairwise intns we solve a system of linear equations
% (Ax=b). The 'A'(Larger_M) and 'b'(The_big_R_matrix) are obtained below

block=zeros(1,n);
Small_M=zeros(n-1);
Large_M=[];
for i=1:n
    Small_M=zeros(n-1,n);
    for j=1:n-1
        Small_M(j,n*(j-1)+1:n*j)=loo_abun(i,:);
    end
    counter=1;
    flag=0;
    Medium_M=zeros(n-1,n*n);
    for k=1:n
        if k==i
            Medium_M(:,(counter-1)*n+1:(counter)*n)=zeros(n-1,n);
            counter=counter+1;
            flag=1;
        else
            Medium_M(:,(counter-1)*n+1:(counter)*n)=Small_M(:,(counter-flag-1)*n+1:(counter-flag)*n);
            counter=counter+1;
        end
    end
    Large_M=[Large_M;Medium_M]; 
end

Last_appendment=zeros(n-1,n*n);
for i=1:n
    Last_appendment(i,n*(i-1)+i)=1;
end

Larger_M=[Large_M;Last_appendment];

Capital_R=[];
for i=1:n
    int_R=r;
    int_R(i)=[];
    Capital_R=[Capital_R;int_R];
end

The_big_R_matrix=[Capital_R;-1*self_intn];

%obtaining effective pairwise interactions 
intn_parameter_vect=(Larger_M)\The_big_R_matrix;

%rearranging effective pairwise interactions into matrix form
for i=1:n
    for j=1:n
        int_par(i,j)=intn_parameter_vect(n*(i-1)+j); %#ok<*SAGROW>
    end
end

%Predicting abundances in the n-member community
abundance_n_member_community=(int_par)\r;


% Printing the outputs
disp('Effective Pairwise Interactions estimated using EPICS');
int_par
disp('Predicted Abundances in the original community using EPICS');
abundance_n_member_community'
