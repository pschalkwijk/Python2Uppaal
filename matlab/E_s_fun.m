%% Calculating the E_s matrix for all regions --> E_s x>=0 \forall x \in R_s
function EE = E_s_fun(n,m)

switch n
    case 2
        %   --> Note: first, we calculate E_s for half space and then shift it to
        %   other half space, because of symmetry that exist between region s and
        %   region s+m (E_s = -E_{s+m})
        q = m;         %half of the number of regions 

        % Calculating all the regions bounds in 2D
        thet_region = (1:m);
        thet_min = (thet_region-ones(size(thet_region)))*pi/m;
        thet_max = thet_region*pi/m;

        EE = cell(2*q,1);   %note: this is 2*q (calculating E_s for all region!)
        for j=1:q

            clear sign_Val a_1 a_2 % Note: a_1^T x>=0 and a_2^T x>=0
            if (thet_min(1,j)<pi/2)  
                sign_Val = 1;
            end

            if (thet_min(1,j)>=pi/2)
                sign_Val = -1;
            end

            a_1 = sign_Val*[cos(thet_min(1,j)) -sin(thet_min(1,j))]';
            a_2 = sign_Val*[-cos(thet_max(1,j)) sin(thet_max(1,j))]';
            EE{j,1} = [a_1';a_2'];
            EE{j+m,1} = -EE{j,1}; %because of symmetry we mentioned before.
        end
        
    otherwise
        
        thet_region = (1:m); %pi/m is theta_i's increment for each region
        sets = cell(1,n-1);
        for i=1:n-1
            sets{1,i} = thet_region; 
        end
        c = cell(1,numel(sets));
        [c{:}] = ndgrid(sets{:});
        permu_All = cell2mat(cellfun(@(v)v(:),c,'UniformOutput',false))';
        permu_sec = [permu_All(1:n-2,:);permu_All(n-1,:)+m];
        permu_second = [permu_All permu_sec];
        
        q = m^(n-1);            % Here, we are using the symmetry property 
        EE = cell(2*q,1);       % but in higher dimensions. But it is a 
                                % little bit more tricky!!
        
        thet_min = (permu_All-ones(size(permu_All)))*pi/m;
        thet_max = (permu_All)*pi/m;
        
        for ww=1:q
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            % computing the middle point of each region
            thet_mid = (thet_min(:,ww)+thet_max(:,ww))/2;
            x_mid = zeros(n,1);
            switch n
                case 3
                    x_mid(1,1) = cos(thet_mid(1,1));
                    x_mid(2,1) = prod(sin(thet_mid));
                    x_mid(3,1) = sin(thet_mid(1,1))*cos(thet_mid(2,1));
                otherwise
                        x_mid(1,1) = cos(thet_mid(1,1));
                        for r=2:n-2
                            x_mid(r,1) = prod(sin(thet_mid(1:r-1,1)))*cos(thet_mid(r,1));
                        end
                        x_mid(n-1,1) = prod(sin(thet_mid));
                        x_mid(n,1) = prod(sin(thet_mid(1:n-2,1)))*cos(thet_mid(n-1,1));
            end
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            
        % % ---------------------------------------------------------------------%%
        % for the case that one of theta_min_i's is 0.
            if (min(thet_min(1:n-2,ww))==0) || (max(thet_max(1:n-2,ww))==pi)  %%%%%%%%%%%%%%Note this changed wrt intersection function!!!!!
             % ---------------------------------------------------------------------%% 
                Face_thet = cell(n-1,2); % the columns in Face_thet are 
                                         % representer of hyperplanes 
                                         % related to thet_{i,min} and 
                                         % thet_{i,max}
                for i=1:n-1

                    thet_min_reduce = thet_min(:,ww);
                    thet_min_reduce(i,:) = [];
                    thet_max_reduce = thet_max(:,ww);
                    thet_max_reduce(i,:) = [];

                    thet_all_reduce = cell(1,n-2);                                     %%%%%%%%%%%%%%%
                    for w=1:n-2
                        thet_all_reduce{1,w} = [thet_min_reduce(w,1) thet_max_reduce(w,1)];
                    end
                    thet_comb_reduce = cartesianProduct(thet_all_reduce)';

                    %the idea for building Face_thet is given in following:
                    %(n-1) is related
                    %to the number of thet_i's and 2^(n-2) is related to
                    %the fact that we keep thet_i constant for its lower 
                    %and upper hyperplanes.
                    Face_thet{i,1} = zeros(n-1,2^(n-2)); 
                    Face_thet{i,1}(i,:) = Face_thet{i,1}(i,:)+thet_min(i,ww); %<-- this is exactly the idea of fixing thet_{i,min}
                    Face_thet{i,2} = zeros(n-1,2^(n-2));
                    Face_thet{i,2}(i,:) = Face_thet{i,2}(i,:)+thet_max(i,ww); %<-- this is exactlt the idea of fixing thet_{i,max}
                    k = 1;
                    for j=1:n-1
                        if j~=i
                            Face_thet{i,1}(j,:) = thet_comb_reduce(k,:);   %the procedure here is based on the way we defined thet_comb_reduce
                            Face_thet{i,2}(j,:) = thet_comb_reduce(k,:);   %hence, we add the combinations that we ignored thet_i in them.
                            k = k+1;
                        end
                    end
                end
             % ---------------------------------------------------------------------%%       
                Face_X = cell(n-1,2); %each column in blocks of Face_X is a
                %representer of the corresponding extreme rays of the 
                %thet_i and also note that each two horizontal bloack are
                %related to the corresponding min and max hyperplanes of
                %thet_i.
                
                for p=1:n-1
                    for t=1:2
                        Face_X{p,t} = zeros(n,2^(n-2));
                        switch n
                            case 3
                                Face_X{p,t}(1,:) = cos(Face_thet{p,t}(1,:));
                                Face_X{p,t}(2,:) = prod(sin(Face_thet{p,t}));
                                Face_X{p,t}(3,:) = sin(Face_thet{p,t}(1,:)).*...
                                    cos(Face_thet{p,t}(2,:));
                            otherwise
                                    Face_X{p,t}(1,:) = cos(Face_thet{p,t}(1,:));
                                    for r=2:n-2
                                        Face_X{p,t}(r,:) = prod(sin(Face_thet{p,t}(1:r-1,:)))...
                                            .*cos(Face_thet{p,t}(r,:));
                                    end
                                    Face_X{p,t}(n-1,:) = prod(sin(Face_thet{p,t}));
                                    Face_X{p,t}(n,:) = prod(sin(Face_thet{p,t}(1:n-2,:)))...
                                        .*cos(Face_thet{p,t}(n-1,:));
                        end
                    end
                end
             % ---------------------------------------------------------------------%%     
            %   this part is added to take into account the numerical precision of
            %   Matlab. The rounding procedure is highly related to the m variable that
            %   user is using. The rounding variable we propose here is for m=100.
                round_Var = 10^4;
                for z=1:n-1
                    for d=1:2
                        Face_X{z,d} = round_Var*Face_X{z,d};
                        Face_X{z,d} = fix(Face_X{z,d});
                        Face_X{z,d} = Face_X{z,d}/round_Var;
                    end
                end

             % ---------------------------------------------------------------------%%
             %Note: here, from all 2^(n-2) vartices on each min and max
             %hyperplanes for thet_i, there should exist (n-1) unique
             %vertices in order to build that hyperplane if that is not the
             %case, as it is in this case, some hyperplanes in cartesian
             %coordinates collapse.
                for ii=1:n-1
                    Face_X_Unique_1 = (unique((Face_X{ii,1})','rows'))';
                    if size(Face_X_Unique_1,2)<(n-1)
                        Face_X{ii,1} = [];
                    else
                        Face_X{ii,1} = [];
                        Face_X{ii,1} = Face_X_Unique_1(:,1:n-1);
                    end
                    Face_X_Unique_2 = (unique((Face_X{ii,2})','rows'))';
                    if size(Face_X_Unique_2,2)<(n-1)
                        Face_X{ii,2} = [];
                    else
                        Face_X{ii,2} = [];
                        Face_X{ii,2} = Face_X_Unique_2(:,1:n-1);
                    end
                end

            % % ---------------------------------------------------------------------%%
                % computing E_s of each cone --> {\forall x\in R_s, E_s x>=0}           
                num_Face = 0; %number of faces, we have for this region.
                clear i j
                for i=1:(n-1)
                    for j=1:2
                        if (isempty(Face_X{i,j})~=1)
                            num_Face = num_Face+1;
                        end
                    end
                end

                EE_reg = zeros(num_Face,n);
                e_crammer = ones(n,1);
                row_counter = 0;
                for i=1:n-1
                    for j=1:2
                        if (isempty(Face_X{i,j})~=1)
                            row_counter = row_counter+1;
                            E_prim = [(Face_X{i,j})';zeros(1,n)];
                            E_prim(:,end) = E_prim(:,end)+e_crammer;%havering plane 
                            normal_vec = zeros(n,1);
                            for k=1:n
                                E_det = E_prim;
                                E_det(:,k) = e_crammer;
                                normal_vec(k,1) = det(E_det);
                            end
                            if (normal_vec'*x_mid<0)
                                normal_vec = -normal_vec;
                            end
                            normal_vec = normal_vec/norm(normal_vec);
                            EE_reg(row_counter,:) = normal_vec';
                        end
                    end
                end
                EE{ww,1} = EE_reg;
                
                % finding the other region with -E_s as its inequality matrix
                rr = zeros(n-1,1);

                for h=1:n-2
                    if permu_second(h,ww)<=(m/2)
                        rr(h,1) = permu_second(h,ww)+...
                            2*(m/2-permu_second(h,ww))+1;
                    else
                        rr(h,1) = permu_second(h,ww)-...
                            2*(permu_second(h,ww)-m/2)+1;
                    end
                end
                rr(n-1,1) = permu_second(n-1,ww)+m;
                %finding the region based on the permu second $$ NOTE $$
                ss = find(ismember(permu_second',rr','rows'));
                
                EE{ss,1} = -EE_reg;  % NOTE: VERY IMPORTANT --> Because for   NOTE --> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                     % each column of permu_second, we are
                                     % finding the E_s, here, it should be
                                     % like this. Later, also, we can check
                                     % it whether it is through or not
                                     
            % % ---------------------------------------------------------------------%%    
            else
            % % ---------------------------------------------------------------------%%

             % ----------------none of thet_i is 0 or pi (i={1:n-2})----------------%%
             % deriving the angle combinations for the rays we want to consider without
             % taking into account the cases A_{cl,d} losses rank.
                Face_thet = cell(n-1,2);%one row for each theta_i containig two blocks -/+ 
                for i=1:n-1
                    Face_thet{i,1} = zeros(n-1,n-1);
                    Face_thet{i,1}(:,1) = thet_min(:,ww);
                    Face_thet{i,2} = Face_thet{i,1};
                    Face_thet{i,2}(i,1) = Face_thet{i,2}(i,1)+pi/m;
                    k = 2; %variable introduced to allocate the vecotrs inside each Face_thet_i
                           %the reason it starts from 2 is because the
                           %first column is already assigned.
                    for j=1:(n-1)
                        if j~=i
                            Face_thet{i,1}(:,k) = Face_thet{i,1}(:,1);
                            Face_thet{i,1}(j,k) = Face_thet{i,1}(j,k)+pi/m;
                            Face_thet{i,2}(:,k) = Face_thet{i,2}(:,1);
                            Face_thet{i,2}(j,k) = Face_thet{i,2}(j,k)+pi/m;
                            k = k+1;
                        end
                    end
                end

            % % ---------------------------------------------------------------------%%
            % deriving the extreme rays for each Face_x without
            % taking into account the cases A_{cl,d} losses rank.
                Face_X = cell(n-1,2);
                for p=1:n-1
                    for t=1:2
                        Face_X{p,t} = zeros(n,n-1);
                        switch n
                            case 3
                                Face_X{p,t}(1,:) = cos(Face_thet{p,t}(1,:));
                                Face_X{p,t}(2,:) = prod(sin(Face_thet{p,t}));
                                Face_X{p,t}(3,:) = sin(Face_thet{p,t}(1,:)).*...
                                    cos(Face_thet{p,t}(2,:));
                            otherwise
                                    Face_X{p,t}(1,:) = cos(Face_thet{p,t}(1,:));
                                    for r=2:n-2
                                        Face_X{p,t}(r,:) = prod(sin(Face_thet{p,t}(1:r-1,:)))...
                                            .*cos(Face_thet{p,t}(r,:));
                                    end
                                    Face_X{p,t}(n-1,:) = prod(sin(Face_thet{p,t}));
                                    Face_X{p,t}(n,:) = prod(sin(Face_thet{p,t}(1:n-2,:)))...
                                        .*cos(Face_thet{p,t}(n-1,:));
                        end
                    end
                end

            % % ---------------------------------------------------------------------%%
                % computing E_s of each cone --> {\forall x\in R_s, E_s x>=0}           
                EE_reg = zeros(2*(n-1),n);  %2*(n-1) because in normal 
                                            % each cone ends up with this
                                            % number of hyperplanes.
                e_crammer = ones(n,1);
                for i=1:n-1
                    for j=1:2
                        E_prim = [(Face_X{i,j})';zeros(1,n)];
                        E_prim(:,end) = E_prim(:,end)+e_crammer;
                        normal_vec = zeros(n,1);
                        for k=1:n
                            E_det = E_prim;
                            E_det(:,k) = e_crammer;
                            normal_vec(k,1) = det(E_det);
                        end
                        if (normal_vec'*x_mid<0)
                            normal_vec = -normal_vec;
                        end
                        normal_vec = normal_vec/norm(normal_vec);
                        EE_reg((i-1)*2+j,:) = normal_vec';
                    end
                end
                EE{ww,1} = EE_reg;
                % finding the other region with -E_s as its inequality matrix
                rr = zeros(n-1,1);

                for h=1:n-2
                    if permu_second(h,ww)<=(m/2)
                        rr(h,1) = permu_second(h,ww)+...
                            2*(m/2-permu_second(h,ww))+1;
                    else
                        rr(h,1) = permu_second(h,ww)-...
                            2*(permu_second(h,ww)-m/2)+1;
                    end
                end
                rr(n-1,1) = permu_second(n-1,ww)+m;
                %finding the region based on the permu second $$ NOTE $$
                ss = find(ismember(permu_second',rr','rows'));
                
                EE{ss,1} = -EE_reg;
            end
        end
end
end















            
            
            
            
            