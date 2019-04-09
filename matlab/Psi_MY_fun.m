% Isontrpic covering %%
% Psi_s: Psi=[Psi_1 ,..., Psi_q] such that x^T Psi_s x >=0 if x in
% R_s ---- del_theta = pi/m;             %region angle increment

function Psi = Psi_MY_fun(n,m)

q = m^(n-1);                  %number of regions 

%%%%%%%%%%%%%%%%%%%% % Calculating all the regions
thet_region = (1:m);
sets = cell(1,n-1);
for i=1:n-1
    sets{1,i} = thet_region; 
end
c = cell(1, numel(sets));
[c{:}] = ndgrid( sets{:} );
% c{1,2}
permu_All = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) )';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thet_min = permu_All-ones(size(permu_All));
thet_min = thet_min*pi/m;
thet_max = permu_All*pi/m;

Psi = cell(q,2*(n-1));

for j=1:q
    clear sign_Val
    switch n
        %%%%%%%%%%%%%%%%%%%%%%%    n=2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 2
            clear sign_Val
            if (thet_min(1,j)<pi/2), sign_Val = 1;end
            if (thet_max(1,j)>pi/2), sign_Val = -1;end

            Psi{j,1} = sign_Val*[(cos(thet_min(1,j)))^2 0;...
                                 0 -(sin(thet_min(1,j)))^2];
            
            Psi{j,2} = sign_Val*[-(cos(thet_max(1,j)))^2 0;...
                                 0 (sin(thet_max(1,j)))^2];

            
        %%%%%%%%%%%%%%%%%%%%%%%    n>2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        otherwise

            for i=1:(n-2)
                clear sign_Val
                if (thet_min(i,j)<pi/2), sign_Val = 1;end
                if (thet_max(i,j)>pi/2), sign_Val = -1;end
                
                Psi{1,2*i-1} = [zeros(i,n);zeros(n-i,i) cos(thet_min(i,j))^2*eye(n-i)];   
                Psi{1,2*i-1}(i,i) = -sin(thet_min(i,j))^2;
                Psi{1,2*i-1} = sign_Val*Psi{1,2*i-1};
                Psi{1,2*i} = [zeros(i,n); zeros(n-i,i) -cos(thet_max(i,j))^2*eye(n-i)];   
                Psi{1,2*i}(i,i) = sin(thet_max(i,j))^2;
                Psi{1,2*i} = sign_Val*Psi{1,2*i};
            end
            clear sign_Val
            if (thet_min(n-1,j)<pi/2), sign_Val = 1;end
            if (thet_max(n-1,j)>pi/2), sign_Val = -1;end
            
            Psi{1,2*(n-2)+1} = zeros(n);
            Psi{1,2*(n-2)+1}(n-1:n,n-1:n) = [cos(thet_min(i,j))^2 0;0 -sin(thet_min(i,j))^2];
            Psi{1,2*(n-2)+1} =  sign_Val*Psi{1,2*(n-2)+1};
            Psi{1,2*(n-1)} = zeros(n);
            Psi{1,2*(n-1)}(n-1:n,n-1:n) = [-cos(thet_max(i,j))^2 0;0 sin(thet_max(i,j))^2];
            Psi{1,2*(n-1)} = sign_Val*Psi{1,2*(n-1)};
            

    end

end 






