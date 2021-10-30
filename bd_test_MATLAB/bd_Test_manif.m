function [boundary_points,boundary_index,distances] = bd_Test_manif(point_cloud,test_index,normal_vectors,eps,r)
%%Given a point cloud and estimated normal vectors, computes the points
%%within e distance of the boundary with high probability. Outputs boundary
%%points, indices of the boundary points, and their estimated distance to
%%the boundary. Uses the 2nd order test as described in the paper.

%INPUT:
%point_cloud
%test_index - indices of points to be tested for the boundary
%normal_vectors - estimated normal vectors, normalized to unit length. Normal vectors can be computed using the function estimated_normal
%eps - thickness of the boundary we seek to identify
%r - neighborhood radius (as in the paper)
% type: type of the test. type=1 is the 1st order test, and type=2 is the 2nd order test.


%OUTPUT:
%boundary_points - boundary points of point_cloud
%boundary_index - indices of boundary points, as a subset of point_cloud
%distances - estimated distances of points tested, according to test_index
    
    normal_vectors=normr(normal_vectors);

    boundary_index=[];
    point_cloud_idx=(1:length(point_cloud)).';
    idx=rangesearch(point_cloud,point_cloud,r);
    distances=zeros(length(point_cloud),1);
    for j=1:length(test_index)
        i=test_index(j);
        ref_point=point_cloud(i,:);
        Test=1;
        
        k=find(point_cloud_idx==i);
        
        n0=normal_vectors(i,:);        

        %
        k_idx=(cell2mat(idx(k))).';
        
        normal_idx=point_cloud_idx(k_idx,:);
        
        bpoint_cloud=point_cloud(k_idx,:);
        ball_normals=normal_vectors(normal_idx,:);
        
        bcent=mean(bpoint_cloud,1); %%barycenter
        cov_mat=cov(bpoint_cloud-bcent); %%covariance matrix - ref_point or barycenter?
        [V,D]=eig(cov_mat); %%eigenvectors of covariance matrix
        eigval=diag(D)./sum(diag(D));
       
        dim=length(eigval);

        
        %find eigenvalues larger than 1/10 of the largest one
        meig_idx=find(eigval>0.1*eigval(dim));
                
        A=V(:,meig_idx); %%matrix of basis, columns vectors
        P=(A/(A.'*A))*(A.'); %projection onto the PCA subspace
        %pdet=det(P);
        %P=P/(det(P)^(1/(m^2)));
        %pdet=det(P);
                
        ball_normals=normr((P*ball_normals.').'); %%projected normals of other points
        n0=normr((P*n0.').'); %%projected normal at x0
        ptdiffmat=repmat(ref_point,length(k_idx),1)-bpoint_cloud;
        ptdiffmat=(P*ptdiffmat.').';
        n0mat=repmat(n0,length(k_idx),1);
        
      
        diff_normals=0.5*(ball_normals-n0mat);
        
        sharp1=sum(n0mat.*ball_normals,2);
        scale_sharp=sharp1>0;
        estimator_sharp=sum(ptdiffmat.*(n0mat+(diff_normals.*scale_sharp)),2);
        dist_est_sharp=max(estimator_sharp);
        
        distances(j)=dist_est_sharp;

        if(size(bpoint_cloud,1)<2)%somewhat arbitrary...
           Test=0;
           fprintf("Insufficient number of points");
        else        
                %Test for boundary point
            threshold=(3/2)*eps;
            if(dist_est_sharp>threshold)
            %if(dist_est_firstorder>threshold)
                Test=0;
            end
        end
        
        
        %Extract indices of boundary points
        if (Test==1)
            boundary_index=[boundary_index,i];
        end

        disp([j,length(test_index)]);
   
    end
    
    if(isempty(boundary_index)==1)
        fprintf("Could not identify boundary points")
        boundary_points=point_cloud(boundary_index,:);
    else
        boundary_points=point_cloud(boundary_index,:);
    end
end