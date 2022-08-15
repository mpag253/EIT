function [K,pK,SpK,i_K,M_B,C,D] = generate_fem_model(parameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[params_m,params_e,params_p] = deal(parameters{:});
[nodes,tris,nn,bdy_indx,bdy_elems] = deal(params_m{:});
[n_elec,in_elec,n_per_elec,z_elec,N] = deal(params_e{:});
[sigma] = deal(params_p{:});

[K,pK,SpK,i_K] = get_fem_K(nodes,tris,sigma,nn);
M_B = get_fem_MB(nn,n_elec,nodes,bdy_indx,in_elec,n_per_elec);
C = get_fem_C(bdy_elems,nn,n_elec,bdy_indx,in_elec,nodes,z_elec,N);
D = get_fem_D(n_elec,bdy_indx,in_elec,nodes,N,z_elec);

end


%% SUB-FUNCTIONS

function [D] = get_fem_D(n_elec,bdy_indx,in_elec,nodes,N,z_elec)
    % 
    el_lengths=zeros(n_elec,1);
    % 
    for ii=1:n_elec
        el_length=0;
        indx=bdy_indx(in_elec(:,ii));
        for jj=1:length(bdy_indx(in_elec(:,ii)))-1
            el_length=el_length+norm(nodes(indx(jj+1),:)-nodes(indx(jj),:));
        end
        el_lengths(ii)=el_length;
    end
    % 
    D=1/z_elec*N'*diag(el_lengths)*N;
end

function [C] = get_fem_C(bdy_elems,nn,n_elec,bdy_indx,in_elec,nodes,z_elec,N)
    % elec_nodes=unique(bdy_el);
    elec_nodes=bdy_elems(:,1);
    
    Chat=sparse(nn,n_elec);
    for ii=1:length(elec_nodes)
        for jj=1:n_elec
            inter=0;
            if sum(elec_nodes(ii)==bdy_indx(in_elec(:,jj)))==1
                
                [rr,cc]=find(elec_nodes(ii)==bdy_elems);
                
                if sum(bdy_elems(rr(1),2-mod(cc(1)+1,2))==bdy_indx(in_elec(:,jj)))
                                        
                    Chat(elec_nodes(ii),jj)=Chat(elec_nodes(ii),jj)+1/2*norm(...
                        nodes(bdy_elems(rr(1),2-mod(cc(1)+1,2)),:)-...
                        nodes(bdy_elems(rr(1),cc(1)),:));
                end
                if sum(bdy_elems(rr(2),2-mod(cc(2)+1,2))==bdy_indx(in_elec(:,jj)))
                    Chat(elec_nodes(ii),jj)=Chat(elec_nodes(ii),jj)+1/2*norm(...
                        nodes(bdy_elems(rr(2),2-mod(cc(2)+1,2)),:)-...
                        nodes(bdy_elems(rr(2),cc(2)),:));
                end
            end
    
        end
    end
    C=-1/z_elec*Chat*N;
end

function [M_B] = get_fem_MB(nn,n_elec,nodes,bdy_indx,in_elec,n_per_elec)
    M_B = sparse([],[],0,nn,nn);
    for ii = 1:n_elec
        [pM_b,~,i_Mb] = PreMassCurved1D(nodes(bdy_indx(in_elec(:,ii)),:),0,0);
        M_b = sparse(i_Mb(:,1),i_Mb(:,2),pM_b*ones(n_per_elec(ii),1),n_per_elec(ii),n_per_elec(ii));
        M_B(bdy_indx(in_elec(:,ii)),bdy_indx(in_elec(:,ii))) = M_b;
    end
end

function [K,pK,SpK,i_K] = get_fem_K(nodes,tris,sigma,nn)
    [pK,SpK,i_K] = Prestiff2D(nodes,tris,1);
    % K = sparse(i_K(:,1),i_K(:,2),pK*sigma_true,nn,nn);
    % K = sparse(i_K(:,1),i_K(:,2),pK*exp(sigma_true),nn,nn);
    K = sparse(i_K(:,1),i_K(:,2),pK*exp(sigma.*ones(nn,1)),nn,nn);
end