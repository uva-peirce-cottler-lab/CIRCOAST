function bw_protect = VesselGen_SelectProtectedSegment(bw_skel,varargin);
% ARGUMENT PARSING
% What is not an Parameter is an optinal arg (param words must fail
% validation for add optional).
p = inputParser;
p.addParameter('FigureVisible','on', @(x) strcmp(x,'off') || strcmp(x,'on'));
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

NOWHERE_TO_GO=0;
MADE_IT_TO_FINISH=0;
img_dim =size(bw_skel);

bw_border = false(img_dim);
bw_border(1,:)=1; bw_border(:,1)=1;
bw_border(end,:)=1; bw_border(:,end)=1;



bw_pad_border = padarray(bw_border,[1 1]);
bw_pad_skel = padarray(bw_skel,[1 1]);

bw_edgepoints = bwmorph(bw_pad_skel,'endpoints') & bw_pad_border;
img_ind =1:prod(img_dim+2);
ep_ind = img_ind(bw_edgepoints);
rand_ind = randperm(numel(ep_ind));


for k=1:numel(ep_ind)
    MADE_IT_TO_FINISH=0;
    NOWHERE_TO_GO=0;
    
    % Radomly select start point on edge of skeleton
    start_ind = ep_ind(rand_ind(k));
    [r, c] = ind2sub(img_dim+2, start_ind);
    
    
    % Calculate trave kernel based on position of start
    bw_pad_ctr = false(img_dim+2);
    bw_pad_ctr(round((img_dim(1)+1)/2),round((img_dim(2)+1)/2))=1;
    ed = bwdist(bw_pad_ctr);
    base_kern = ed(r-1:r+1,c-1:c+1);
    trav_kern = max(base_kern(:))-base_kern+1;
    trav_kern=trav_kern./max(trav_kern(:));
    trav_kern(2,2)=0;
    
    % Kernels to keep track of movement
    c_kern = repmat([-1 0 1],[3 1]);
    r_kern = repmat([-1; 0; 1],[1 3]);
    
    
    curr_rc = [r c];
    bw_pad_history =false(img_dim+2);
    bw_pad_history(start_ind)=1;
    % hw = waitbar(0,'Finding Central Vessel...');
    for n=2:10*prod(img_dim+2)
        %Get 9 neighborhood in image
        ri = curr_rc(n-1,1);
        ci = curr_rc(n-1,2);
        
        
        nhood = bw_pad_skel(ri-1:ri+1,ci-1:ci+1) & ...
            ~bw_pad_history(ri-1:ri+1,ci-1:ci+1);
        
        % Multiply by kern figure out which is max
        score = nhood.*trav_kern;
        
        % Find best pixel to move to
        [~,ind]= max(score(:));
        
        % Calculate next move
        curr_rc(n,1:2) = curr_rc(n-1,1:2) + [r_kern(ind(1)) c_kern(ind(1))] ;
        
        % Add move to bw_pad_history
        next_ind = sub2ind(img_dim+2,curr_rc(n,1),curr_rc(n,2));
        bw_pad_history(next_ind)=1;
        
        % Exit Conditions
        if bw_pad_border(next_ind) && n>(img_dim(1)/2)
            MADE_IT_TO_FINISH=1;
            fprintf('Kernel traveled across image.\n')
        end
        
        if all(score==0); 
            NOWHERE_TO_GO=1; 
            fprintf('Kernel Caught in Dead End.\n')
        end
        
        
        if strcmp(FigureVisible,'on') && (rem(n,25)==0 || MADE_IT_TO_FINISH || NOWHERE_TO_GO);
            dil_hist = imdilate(bw_pad_history,strel('disk',0,0));
            dil_skel = imdilate(bw_pad_skel,strel('disk',0,0));
            dil_curr = false(img_dim+2);
            dil_curr(next_ind)=1;
            
            subplot(2,2,1); imshow(dil_hist);
            subplot(2,2,2); imshow(im2uint8(cat(3,dil_hist,dil_curr,dil_skel)));
            subplot(2,2,3); imshow(trav_kern)
            subplot(2,2,4); plot(curr_rc(:,1),curr_rc(:,2))
            pause(.1);
        end
        
        if NOWHERE_TO_GO; break; end
        if MADE_IT_TO_FINISH; break; end
    end
    close(gcf)
    
    if MADE_IT_TO_FINISH; break; end
    
end
fprintf('NOWHERE_TO_GO:%i\n',NOWHERE_TO_GO);
fprintf('MADE_IT_TO_FINISH:%i\n',MADE_IT_TO_FINISH)

 
% if ~any(curr_rc(n,1:2) - 512>=0); keyboard; end
if sum(sum(bw_pad_border & bw_pad_history))<2;keyboard; end

bw_protect = bw_pad_history(2:end-1,2:end-1);
