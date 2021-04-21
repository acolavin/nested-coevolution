function [r2, sectors] = MirrorTree(msasort, secs, filepathin)
%% Adding random sectors

for i = 1:length(secs)
    for j = 1:length(secs{i})
        seclen(j) = length(secs{i}{j});
    end
    seclen(seclen < 3) = NaN;
    
%     sectorSize = [min(seclen) max(seclen)];
%     sectorNum = 20;
%     randType = 1;
    
    sectorSize = seclen;
    sectorNum = 20;
    randType = 3;
    
    randSectors = randomSectors(size(msasort{i}, 2), sectorSize, sectorNum, randType);
    
    sectors{i} = [secs{i} randSectors];
end

r2 = NaN;

%% Calculate distance matrices for all sectors

for i = 1:length(sectors{1})
    if ~exist([filepathin num2str(i)], 'file')
        disp(['Working on sector: ' num2str(i) '/120'])
        tic
        sectormsa1 = msasort{1}(:,sectors{1}{i});
        if ~isempty(sectors{1}{i})
            distmat = calcJukesCantorDist(sectormsa1);        
        else
            num = size(msasort{1}, 1);
            distmat = nan(num,num);
        end
        toc
        save([filepathin num2str(i)], 'distmat', '-v7.3')
        toc
    end
end   

% % loading data
% disp('Making cell of distmat')
% tic
% for j = 1:length(sectors{1})
%     load([filepathin num2str(j)])
%     distmat1{j} = distmat;   
%     toc
% end
    
%     disp('Saving Distance Mat')
%     tic
%     save([filepathin, 'distmat1', '-v7.3')
%     toc
    
% disp('Done calculating distance matrices!')

%% r-squared values

disp('Finding R Values')


for i = 1:length(sectors{1})
    load([filepathin num2str(i)])
    distance_matrix_1 = reshape(triu(distmat), [], 1);
    for j = 1:length(sectors{1})
        disp(['Row: ' num2str(i) '/120'])
        disp(['Col: ' num2str(j) '/120'])
        tic
        
        load([filepathin num2str(j)])
        distance_matrix_2 = reshape(triu(distmat), [], 1);
        
        if isempty(distance_matrix_1) == 1 || isempty(distance_matrix_2) == 1
            r2(i,j) = NaN;
        else
            r2mat = corrcoef(distance_matrix_1, distance_matrix_2);
            r2(i,j) = r2mat(1,2);
        end
        toc
    end
end

%% test figure

% hf = figure;
% imagesc(r2)
% title('r-squared')
% colorbar
% set(gca, 'fontsize', 16)
% set(gcf, 'color', 'w')
% set(gcf, 'Position', [100, 100, 1500, 1000])


%% DONE
% disp('Done with MirrorTree!')

end

function sectorsRand = randomSectors(len, sectorSize, sectorNum, randType)

    if randType == 1
        for i = 1:sectorNum
            sectorsRand{i} = randi([1 len], randi(sectorSize, 1), 1)';
        end
    end
    
    if randType == 2
        for i = 1:sectorNum
            f = [0:randi(sectorSize,1)]+randi([1 len], 1, 1);
            f(f > len) = f(find(f>len,1):end)-len;

            sectorsRand{i} = f;
        end
    end
    
    if randType == 3
        num = 1;
        for i = 1:sectorNum
            for j = 1:5 % number of bootstraps
                
                if ~isnan(sectorSize(i))               
                    sectorsRand{num} = randi([1 len], sectorSize(i), 1)';                
                else
                    sectorsRand{num} = [];                  
                end
                num = num + 1;
                
            end
        end
    end
    
end

% Calculate pairwise Jukes Cantor Distance of a given MSA
function pairwised = calcJukesCantorDist(msa)
    [rows, columns] = size(msa);
    pairwised = zeros(1,rows);
    parfor i = 1:rows
        rowRepeat = repmat(msa(i,:),rows,1);
        matches = rowRepeat ~= msa;
        p = sum(matches,2)/columns;
        p = min(p,.949);
        pairwised(i,:) = -19/20*log(1-p*20/19);
    end
    
end

