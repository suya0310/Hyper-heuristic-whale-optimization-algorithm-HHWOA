function OBLPositions=OBLinitialization(SearchAgents_no,dim,lb,ub,Positions)
OBLPositions=zeros(SearchAgents_no,dim);
Boundary_no=size(ub,2); % numnber of boundaries    
if Boundary_no==1
         for j=1:size(Positions,2)
           OBLPositions(:,j)=lb+ub-Positions(:,j);     
         end  
end
if Boundary_no>1
         for j=1:size(Positions,2)
           OBLPositions(:,j)=lb(j)+ub(j)-Positions(:,j);     
         end
end