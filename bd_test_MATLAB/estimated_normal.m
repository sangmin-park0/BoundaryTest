function nvec = estimated_normal(point_cloud,r)
    %%Compute estimated normal vector using graph laplacian, with weights
    %%normalizing by the degree. The code assumes the smallest array
    %%dimension of point cloud is the dimension of the Euclidean space in which the
    %%points lie.
    
    %Input: point_cloud, r (neighborhood radius)
    
    %Output: nvec (estimated normal vectors)
    
    dim=min(size(point_cloud));
    nvec=zeros(length(point_cloud),dim);
    idx=rangesearch(point_cloud,point_cloud,r);


    idx_small=rangesearch(point_cloud,point_cloud,r/2);
    deg_vec = cellfun(@length, idx_small);
    
    for j=1:length(point_cloud)
        ref_point=point_cloud(j,:);
        
        
        ball_idx=(cell2mat(idx(j))).';
        ball_points=point_cloud(ball_idx,:);
                 
        deg_vec_j=deg_vec(ball_idx);
        weights=1./(deg_vec_j); %%simple normalization by degree
               
        %calculate the estimated normal

        diff_mat=ball_points-repmat(ref_point,length(ball_points(:,1)),1);
        est_normal=sum(weights.*diff_mat,1);
        
        %disp([j,length(point_cloud)]);
        nvec(j,:)=est_normal;
        
    end

end