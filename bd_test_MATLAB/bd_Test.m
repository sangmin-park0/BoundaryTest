function [boundary_points,boundary_index,distances] = bd_Test(point_cloud,test_index,normal_vectors,eps,r,type)
%%Given a point cloud and estimated normal vectors, computes the points
%%within e distance of the boundary with high probability. Outputs boundary
%%points, indices of the boundary points, and their estimated distance to
%%the boundary

%INPUT:
%point_cloud
%test_index - indices of points to be tested for the boundary
%normal_vectors - estimated normal vectors. Normal vectors can be computed using the function estimated_normal
%eps - thickness of the boundary we seek to identify
%r - neighborhood radius (as in the paper)
%type: type of the test. type=1 is the 1st order test, and type=2 is the 2nd order test.


%OUTPUT:
%boundary_points - boundary points of point_cloud
%boundary_index - indices of boundary points, as a subset of point_cloud
%distances - estimated distances of points tested, according to test_index
    
    %normalize the inward normal vectors
    normal_vectors=normr(normal_vectors);
    
    boundary_index=[];
    all_idx=(1:length(point_cloud)).';
    idx=rangesearch(point_cloud,point_cloud,r);
    distances=zeros(length(test_index),1);
    for j=1:length(test_index)
        i=test_index(j);
        ref_point=point_cloud(i,:);
        Test=1;

        
        n0=normal_vectors(i,:);
        

        %
        k_idx=(cell2mat(idx(i))).';
        %all_idx=(1:length(all_points)).';
        %ball_idx=all_idx_input(j_idx,:);
        %ball_points=data(ball_idx,:);
        %ball_normals=all_normals(ball_idx,:);
        
        normal_idx=all_idx(k_idx,:);
        
        ball_points=point_cloud(k_idx,:);
        ball_normals=normal_vectors(normal_idx,:);
        
        ptdiffmat=repmat(ref_point,length(k_idx),1)-ball_points;
        n0mat=repmat(n0,length(k_idx),1);
        
        avg_mat=0.5*(n0mat+ball_normals);

        term1=sum(ptdiffmat.*avg_mat,2);
    

        
       % estimator2=sum(ptdiffmat.*normr(avg_mat),2);
       % dist_est2=max(estimator2);
        
        diff_normals=0.5*(ball_normals-n0mat);
        norm_diff=vecnorm((ball_normals-n0mat).').';
       % square=1-norm_diff.*norm_diff;
       % square=[square.';zeros(1,length(square))].';
       % scale=max(square.').';
        
       % quad=[(1-0.1*norm_diff.^10).';zeros(1,length(square))].';
       % scale_quad=max(quad.').';
       % estimator_quad=sum(ptdiffmat.*(n0mat+diff_normals.*scale_quad),2);
       % dist_est_quad=max(estimator_quad);
        
        
       % estimator_new=sum(ptdiffmat.*(n0mat+diff_normals.*scale),2);
       % dist_est_new=max(estimator_new);
        
        sharp1=sum(n0mat.*ball_normals,2);
        scale_sharp=sharp1>0;
        estimator_sharp=sum(ptdiffmat.*(n0mat+diff_normals.*scale_sharp),2);
        dist_est_sharp=max(estimator_sharp);
        
        
        
        switch type
            case 1
                estimator_firstorder=sum(ptdiffmat.*n0mat,2);
                dist_est_firstorder=max(estimator_firstorder);
                distances(j)=dist_est_firstorder;
                dist_est=dist_est_firstorder;
            case 2
                distances(j)=dist_est_sharp;
                dist_est=dist_est_sharp;
        end
       % estimator_firstorder=sum(ptdiffmat.*n0mat,2);
       % dist_est_firstorder=max(estimator_firstorder);
        
        
        %dist_true=pdist2(ref_point,[1,1,1])-0.5; %invsphere
        
        if(size(ball_points,1)<5)%somewhat arbitrary...
           Test=0;
           fprintf("Insufficient number of points");
        else        
                %Test for boundary point
            threshold=(3/2)*eps;
            if(dist_est>threshold)
               Test=0;
            end
        end
        
        
        %Extract indices of boundary points
        if (Test==1)
            boundary_index=[boundary_index,i];
        end
        %disp(boundary_index);
        disp([j,length(test_index)]);
   
    end
    
    if(isempty(boundary_index)==1)
        fprintf("Could not identify boundary points")
        boundary_points=point_cloud(boundary_index,:);
    else
        boundary_points=point_cloud(boundary_index,:);
    end
end