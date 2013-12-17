function hist = getChainCodeHist(images)
    hist = [length(images), 8];
    for i = 1:length(images)
        hist(i, 1:8) = getChainCodeHistOfOne(images{i})';
    end
end

function hist = getChainCodeHistOfOne(image)
    
    [x,y] = find(image == 1, 1);
    image(x, y) = 0;
    hist = zeros([8, 1]);
    while(sum(sum(image)) ~= 0)
        for k=1:size(image)
            % Top left
            if (image(x - 1, y - 1) == 1)
                hist(1) = hist(1) + 1;
                x = x - 1;
                y = y - 1;
            % Top center
            elseif (image(x, y - 1) == 1)
                hist(2) = hist(2) + 1;
                y = y - 1;
            % Top right
            elseif (image(x + 1, y - 1) == 1)
                hist(3) = hist(3) + 1;
                x = x + 1;
                y = y - 1;
            % Right center
            elseif (image(x + 1, y) == 1)
                hist(4) = hist(4) + 1;
                x = x + 1;
            % Bottom right
            elseif (image(x + 1, y + 1) == 1)
                hist(5) = hist(5) + 1;
                x = x + 1;
                y = y + 1;
            % Bottom center
            elseif (image(x, y + 1) == 1)
                hist(6) = hist(6) + 1;
                y = y + 1;
            % Bottom left
            elseif (image(x - 1, y + 1) == 1)
                hist(7) = hist(7) + 1;
                x = x - 1;
                y = y + 1;
            % Left center
            elseif (image(x - 1, y) == 1)
                hist(8) = hist(8) + 1;
                x = x - 1;
            else
                [x,y] = find(image == 1, 1);
            end
            image(x, y) = 0;
        end
        
    end

%     hist = zeros([8, 1]);
%     
%    for j=1:size(edges)
%         x = edges(j, 1);
%         y = edges(j, 2);
% 
%         for k=1:size(edges)
%             % Top left
%             if (j ~= k && x == edges(k, 1) - 1 && y == edges(k, 2) - 1)
%                 hist(1) = hist(1) + 1;
%             end
% 
%             % Top center
%             if (j ~= k && x == edges(k, 1) && y == edges(k, 2) - 1)
%                 hist(2) = hist(2) + 1;
%             end
% 
%             % Top right
%             if (j ~= k && x == edges(k, 1) + 1 && y == edges(k, 2) - 1)
%                 hist(3) = hist(3) + 1;
%             end
% 
%             % Right center
%             if (j ~= k && x == edges(k, 1) + 1 && y == edges(k, 2))
%                 hist(4) = hist(4) + 1;
%             end
% 
%             % Bottom right
%             if (j ~= k && x == edges(k, 1) + 1 && y == edges(k, 2) + 1)
%                 hist(5) = hist(5) + 1;
%             end
% 
%             % Bottom center
%             if (j ~= k && x == edges(k, 1) && y == edges(k, 2) + 1)
%                 hist(6) = hist(6) + 1;
%             end
% 
%             % Bottom left
%             if (j ~= k && x == edges(k, 1) - 1 && y == edges(k, 2) + 1)
%                 hist(7) = hist(7) + 1;
%             end
% 
%             % Left center
%             if (j ~= k && x == edges(k, 1) - 1 && y == edges(k, 2))
%                 hist(8) = hist(8) + 1;
%             end
% 
%         end
%    end
end