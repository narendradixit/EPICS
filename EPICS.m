% EPICS=Effective Pairwise Interactions for predicting Community Structure

%The following function takes n(number of species in the community),
%species abundances in monoculture and leave-one-out subcommunities as inputs and
%produces effective pairwise interactions and community structure using EPICS as
%outputs. In this code, we assumed that the microbial community follows
%the generalized Lotka-Volterra model. Other population dynamics models can be
%used instead.

% To obtain effective pairwise intns we solve a system of linear equations
% (Ax=b). The 'A'(Larger_M) and 'b'(The_big_R_matrix) are obtained below.
% The algorithm solves n^2 equations of which n*(n-1) equations correspond
% to species abundances in leave-one-out communities and n equations to
% abundances in monocultures. The matrix 'A' is constructed in blocks. The
% first n blocks are of n-1 rows each, with each block corresponding to a
% leave-one-out experiment, starting from species 1 and ending with species n. Within a
% block, each row corresponds to the abundance of one species present in that leave-one-out experiment.
% The last block contains n equations, with each corresponding to a
% monoculture. The unknowns, which are the effective interaction
% coefficients, are arranged in ascending order of species number (a11, a12, etc.). 
% The right-hand side vector 'b' is derived from the monoculture intrinsic
% growth rates.

clear all;
close all;
clc;
disp('Welcome to EPICS...');

%number of species in the community
number_of_species=input('\nInput number of species in the community=');


%monoculture abundances
ind_abun=input('\nInput monoculture abundances \n(Put your entries wrapped in square brackets. \nPut a space or a comma after each entry)=');

%leave-one-out abundances
loo_abun=input('\nInput abundances in leave-one-out subcommunities \n(Put your entries wrapped in square brackets. \n(Rows are subcommunities, columns are species)=');

global n %number of species in the community


n=number_of_species;

%intrinsic growth rate vector
r=1*ones(n,1);

for i=1:n
diagonal(i,1)=1/ind_abun(1,i);%the diagonal elements or self interaction terms are determined from monoculture abundances
end


block=zeros(1,n);
Small_M=zeros(n-1);
Large_M=[];

%this creates matrix 'A'
for i=1:n
    %to create Medium_M, we first create Small_M 
    
    Small_M=zeros(n-1,n);
    for j=1:n-1
        Small_M(j,n*(j-1)+1:n*j)=loo_abun(i,:);
    end
    counter=1;
    flag=0;
    Medium_M=zeros(n-1,n*n);
    
    %this loop combines Small_M's to create Medium_M by adding zeros to
    %account for the missing species across leave-one-out cultures
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


%this creates vector 'b'
Capital_R=[];
for i=1:n
    int_R=r;
    int_R(i)=[];
    Capital_R=[Capital_R;int_R];
end

The_big_R_matrix=-1*[Capital_R;diagonal];

%this computes effective pairwise interactions 
intn_parameter_vect=(Larger_M)\The_big_R_matrix;

%this rearranges effective pairwise interactions into matrix form
for i=1:n
    for j=1:n
        int_par(i,j)=intn_parameter_vect(n*(i-1)+j);
    end
end

%int_par=-1*int_par;

%this predicts abundances in the n-member community
abundance_n_member_community=-1*(int_par)\r;


%this prints the outputs
disp('Effective Pairwise Interactions estimated using EPICS');
int_par
disp('Predicted Abundances in the original community using EPICS');
abundance_n_member_community'
