%% Region determination of a state
% Note: if a state lives on the boundaries of two intersection then the
% output vector y will have 2^(n-1) positive entries.

function s = reg_Det(x,m,permu_All)
% x := state that we want to determine its associated region
n = length(x');
r = norm(x);
switch n
    case 2
        b = 0;      %reminder variable to fix the variable s later
        if x(1,1)<0
            x = -x;
            b = 1;  %means x(1,1)<0
        end
        theta = acos(x(2,1)/r);
        s = ceil(theta*m/pi);
        if s==0
            s = 1;
        end
        if b==1
            s = s+m;
        end
    %case 3
%         b = 0;      %reminder variable to fix the variable s later
%         if x(2,1)<0
%                 x = -x;
%                 b = 1;
%         end
%         theta = zeros(n-1,1);
%         s_reg = zeros(n-1,1);       
%         theta(1,1) = acos(x(1,1)/r);
%         theta(2,1) = acos(x(3,1)/(r*sin(theta(1,1)))); % Why ???
%         s_reg(1,1) = ceil(theta(1,1)*m/pi);
%         if s_reg(1,1)==0
%             s_reg(1,1) = 1;
%         end
%         s_reg(2,1) = ceil(theta(2,1)*m/pi);
%         if s_reg(2,1)==0
%             s_reg(2,1) = 1;
%         end
%         if b==1  %Note: all theta_i's will be affected.
%             if s_reg(1,1)<=(m/2)
%                 s_reg(1,1) = s_reg(1,1)+2*(m/2-s_reg(1,1))+1;
%             else
%                 s_reg(1,1) = s_reg(1,1)-2*(s_reg(1,1)-m/2)+1;
%             end
%             s_reg(2,1) = s_reg(2,1)+m;
%         end
%         s = find(ismember(permu_All',s_reg,'rows')); %finding the region % CHOP: changed s_reg' to s_reg
%     otherwise
%         b = 0;      %reminder variable to fix the variable s later
%         if x(n-1,1)<0
%                 x = -x;
%                 b = 1;
%         end
%         theta = zeros(n-1,1);
%         s_reg = zeros(n-1,1);
%         theta(1,1) = acos(x(1,1)/r);
%         for j=2:n-2
%             theta(j,1) = acos(x(j,1)/(r*prod(sin(theta(1:j-1,1)))));
%         end
%         theta(n-1,1) = acos(x(n,1)/(r*prod(sin(theta(1:n-2,1)))));       
%         
%         for e=1:(n-1)
%             s_reg(e,1) = ceil(theta(e,1)*m/pi);
%             if s_reg(e,1)==0
%                 s_reg(e,1) = 1;
%             end
%         end
%         if b==1
%             for bb=1:n-2
%                 if s_reg(bb,1)<=(m/2)
%                     s_reg(bb,1) = s_reg(bb,1)+2*(m/2-s_reg(bb,1))+1;
%                 else
%                     s_reg(bb,1) = s_reg(bb,1)+2*(s_reg(bb,1)-m/2)+1;
%                 end
%             end
%             s_reg(n-1,1) = s_reg(n-1,1)+m;
%         end
%         s = find(ismember(permu_All',s_reg','rows')); %finding the region
        
        % ALTERNATIVE WAY TO DETERMINE REGION (CHOP EDIT): Use E_s_fun !!!
    otherwise
        E = E_s_fun(n,m);
        for i = 1:2*m^(n-1)
            if E{i,1}*x >= 0
                s = i;
                break
            end
        end
       
end

end
    



     