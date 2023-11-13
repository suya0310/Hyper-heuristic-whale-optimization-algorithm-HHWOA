
% Main paper: % The Hyper-heuristic Whale Optimization Algorithm (HHWOA)
% A Hybrid Hyper-Heuristic Whale Optimization Algorithm for Reusable Launch Vehicle Reentry Trajectory Optimization.
% Su, Ya, Dai, Ying, and Liu, Yi
% Aerospace Science and Technology, Vol. 119, 2021, p. 107200.
% https://doi.org/10.1016/j.ast.2021.107200.

% I acknowledge that this version of WOA has been written using
% a large portion of the following code:
% The Whale Optimization Algorithm.
% Mirjalili, S., and Lewis, A
% Advances in Engineering Software, Vol. 95, 2016, pp. 51–67. 
% https://doi.org/10.1016/j.advengsoft.2016.01.008.

% Harris hawks optimization: Algorithm and applications
% Ali Asghar Heidari, Seyedali Mirjalili, Hossam Faris, Ibrahim Aljarah, Majdi Mafarja, Huiling Chen
% Future Generation Computer Systems, 
% DOI: https://doi.org/10.1016/j.future.2019.02.028
% ___________________________________________________

function [Leader_score,Leader_pos,Convergence_curve,Total_time]=HHWOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
disp('HHWOA is now tackling your problem')
tic
Leader_pos=zeros(1,dim);
Leader_score=inf; 
Positions=initialization(SearchAgents_no,dim,ub,lb);
WOA_Positions=Positions;
Convergence_curve=zeros(1,Max_iter); 
Tfitness=zeros(2*SearchAgents_no,1);

for i=1:SearchAgents_no               
          Ifitness=fobj(Positions(i,:));   
       if Ifitness<Leader_score % Change this to > for maximization problem  
          Leader_score=Ifitness;      % Update alpha 
          Leader_pos=Positions(i,:);
       end         
end 
w=3;
p=rand;
t=0;

while t<Max_iter       
    
     p=abs(cos(w.*acos(p))); 
     OBLPositions=OBLinitialization(SearchAgents_no,dim,lb,ub,Positions);      %  obl    
   % Return back the search agents that go beyond the boundaries of the search space    
   for i=1:size(OBLPositions,1) 
      Flag4ub=OBLPositions(i,:)>ub;
      Flag4lb=OBLPositions(i,:)<lb;
      OBLPositions(i,:)=(OBLPositions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
   end
    
    TPositions=zeros(2*SearchAgents_no,dim);
    for i=1:SearchAgents_no        
       TPositions(i,:)=Positions(i,:);   
       TPositions(i+SearchAgents_no,:)=OBLPositions(i,:);         
    end
   % the fitness of the total solutions 
    for i=1:2*SearchAgents_no        
      Tfitness(i,1)=fobj(TPositions(i,:));          
    end 

    [~,Tindex]=sort(Tfitness);       %求极小值 从小到大排序      
    for newindex=1:SearchAgents_no
       Positions(newindex,:)=TPositions(Tindex(newindex),:);       
    end     
   [Positions,Leader_score,Leader_pos]=DE(Positions,Leader_score,Leader_pos,lb,ub,fobj); 
    
%   When solving an optimal control problem, the following code of smoothing technique can be used
%     for i=1:SearchAgents_no
%          for j=2:dim-1
%              if ((Positions(i,j)-Positions(i,j-1))*(Positions(i,j+1)-Positions(i,j)))<0
%                 Positions(i,j)=0.5*(Positions(i,j-1)+Positions(i,j+1));
%              end
%          end
%     end  
    
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper     
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)        
%         p = rand();        % p in Eq. (2.6)        
        for j=1:size(Positions,2)            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    WOA_Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    WOA_Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end                
            elseif p>=0.5              
               Distance2Leader=abs(Leader_pos(j)-Positions(i,j));               
               WOA_Positions(i,j)=Distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);                
            end            
        end
        % Return back the search agents that go beyond the boundaries of the search space    
        Flag4ub=WOA_Positions(i,:)>ub;
        Flag4lb=WOA_Positions(i,:)<lb;
        WOA_Positions(i,:)=(WOA_Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        if fobj(WOA_Positions(i,:))<fobj(Positions(i,:))
                   Positions(i,:)=WOA_Positions(i,:);    
                   if fobj(Positions(i,:))<Leader_score
                        Leader_score=fobj(Positions(i,:));
                        Leader_pos=Positions(i,:);
                   end 
        end          
    end
    t=t+1;
    Convergence_curve(t)=Leader_score;    
    display(['Iteration:',num2str(t), '   Leader_score:',num2str(Leader_score)]);   
    
end
Total_time=toc;
end

function [Positions,Leader_score,Leader_pos]=DE(Positions,Leader_score,Leader_pos,lb,ub,fobj)
F=0.5;
CR=0.9; 
DEMax_iter=5;
DEPositions=zeros(size(Positions));
t=0;% Loop counter
while t<DEMax_iter
    for i=1:size(Positions,1)           
            kkk=randperm(size(Positions,1));
            kkk(i==kkk)=[];
            jrand=randi(size(Positions,2));  
            for j=1:size(Positions,2)        
                  if (rand<=CR)||(jrand==j)
                   DEPositions(i,j)=Positions(kkk(1),j)+F*(Positions(kkk(2),j)-Positions(kkk(3),j));
%                  DEPositions(i,j)=Leader_pos(j)+F*(Positions(kkk(1),j)-Positions(kkk(2),j));
%                  DEPositions(i,j)=Leader_pos(j)+F*(Positions(kkk(1),j)-Positions(i,j))+F*(Positions(kkk(2),j)-Positions(kkk(3),j));
                  else
                   DEPositions(i,j)=Positions(i,j);             
                  end
            end           
            % Return back the search agents that go beyond the boundaries of the search space       
             Flag4ub=DEPositions(i,:)>ub; Flag4lb=DEPositions(i,:)<lb;             
             DEPositions(i,:)=(DEPositions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;        
             if fobj(DEPositions(i,:))<fobj(Positions(i,:))
                   Positions(i,:)=DEPositions(i,:);    
                   if fobj(Positions(i,:))<Leader_score
                        Leader_score=fobj(Positions(i,:));
                        Leader_pos=Positions(i,:);
                   end 
             end               
    end     
t=t+1; 
end
end





