function bw_vessels = GenerateVessels(img_size, num_vessels, vessel_diam)



img_size = [512 512];
num_vessels = 100;
vessel_diam = 0;


bw_border = false(img_size);
border_ind = [...
    sub2ind(img_size,(1:img_size(1))',ones(img_size(1),1)) ...
    sub2ind(img_size,(1:img_size(1))',img_size(1)*ones(img_size(1),1)) ... 
    sub2ind(img_size,ones(img_size(2),1),(1:img_size(2))') ...
    sub2ind(img_size,img_size(2)*ones(img_size(2),1),(1:img_size(2))')];


bw_border(1,:)=1; bw_border(:,1)=1; bw_border(end,:)=1; bw_border(:,end)=1;
border_ind = find(bw_border);
% img_ind = reshape(1:prod(img_size),img_size);
% bonder_ind = img_ind(bw_border);


% Generate starting coordinates for all vessels along border of image
border_inds = randi(numel(border_ind), [1 num_vessels]);
img_inds = border_ind(border_inds);
[r0, c0] = ind2sub(img_size, img_inds);

% Plot coordinates to confirm they originate on all borders
% plot(r0,c0,'rx');

% Maps movement kernel to degree of direction
kernel_row = [1 1 1; 2 2 2; 3 3 3] - 2;
kernel_col = [1 2 3; 1 2 3; 1 2 3] - 2;

% Temp image to generate current vessel, add to bw_vessels each iteration
temp = false(img_size);
bw_vessels = false(img_size);

for n = 1 : num_vessels
    ri = r0(n); % 505
    ci = c0(n); % 1
    fprintf('Vessel %.0f: %.0f %.0f\n', n, ri, ci);
    
    % Place a kernel where this spawn point is, calculated the row and
    % column coordinates of this kernel
    coord_row_kernel = kernel_row + ri;
    coord_col_kernel = kernel_col + ci;
    
    dist_kernel = sqrt((coord_col_kernel - img_size(2)/2)^2 + ...
        (coord_row_kernel - img_size(2)/2)^2);
    norm_dist_kernel =  1  - (dist_kernel - min(dist_kernel(:))) / ...
        max(max(dist_kernel - min(dist_kernel(:))));
    
    % Rank norm_dist_kernel highest to lowest
    [~, ib] = sort(norm_dist_kernel(:), 'descend');
    [~, ia] = sort(ib);
    iv = 1:numel(norm_dist_kernel);
    sort_kernel = reshape(iv(ia), size(dist_kernel));
    
    % if blah blah set to <= 4
    move_kernel = (sort_kernel <= 5) .* norm_dist_kernel;
    move_kernel(2,2) = 0;
    
%     fcn_params = rand()
    %
    
    rv = ri; cv = ci; 
    for k = 1:(512*5)
        [rshift, cshift] = VesselGenerator_CalculateMove(move_kernel>0);
        rv(k+1) = rv(k) + rshift;
        cv(k+1) = cv(k) + cshift;
        
        if rv(end) == 0 || cv(end) == 0 || ...
                rv(end) == (img_size(1) + 1) || cv(end) == (img_size(2) + 1);
% %            keyboard
%             if ri==1 || ri ==img_size(1);keyboard; end
%             if numel(rv) < 300; keyboard; end
            break;
        end
    end
    rv(end) = []; cv(end) = [];
    temp(sub2ind(img_size, rv,cv)) = 1;
    temp_dil = imdilate(temp,strel('disk',ceil(vessel_diam/2),0));
    
    bw_vessels = bw_vessels | temp_dil;
    temp(:)=0;
end
imshow(bw_vessels)
keyboard
