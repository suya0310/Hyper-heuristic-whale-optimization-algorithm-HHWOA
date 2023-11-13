% Main paper: % The Hyper-heuristic Whale Optimization Algorithm (HHWOA)
% A Hybrid Hyper-Heuristic Whale Optimization Algorithm for Reusable Launch Vehicle Reentry Trajectory Optimization.
% Su, Ya, Dai, Ying, and Liu, Yi
% Aerospace Science and Technology, Vol. 119, 2021, p. 107200.
% https://doi.org/10.1016/j.ast.2021.107200.

close all;clear; clc
rng('shuffle')

Function_name='F10';                % Name of the test function 
SearchAgents_no=30;                 % Number of search agents  30
Max_iteration=500;                  % Maximum numbef of iterations  500
% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
[Rabbit_Energy,Rabbit_Location,CNVG]=HHWOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

display(['The best location is: ', num2str(Rabbit_Location)]);
display(['The best fitness is: ', num2str(Rabbit_Energy)]);

%Draw objective space
figure,
hold on
semilogy(CNVG,'Color','b','LineWidth',1);
title('Convergence curve','FontName','Times New Roman','FontSize',13)
xlabel('Iteration','FontName','Times New Roman','FontSize',16);
ylabel('Best fitness obtained so far','FontName','Times New Roman','FontSize',16);
legend ('HHWOA','FontName','Times New Roman','FontSize',13)
set (gca,'position',[0.11,0.16,0.83,0.75] );   %0.81   0.78
set(gcf,'unit','centimeters','position',[10,10,15,10])  %17  11
set(gca, 'LineWidth',1.0)
set(gca,'FontSize',15,'FontName','Times New Roman');
% axis tight
grid on
axis tight
grid on
box on

        



