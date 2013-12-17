function hist = getChainCodeHist(edges)
    hist = [length(edges), 8];
    for i = 1:length(edges)
        hist(i, 1:8) = getChainCodeHistOfOne(edges{i})';
        disp('Image: ' + i);
    end
end

function hist = getChainCodeHistOfOne(edges)
    hist = zeros([8, 1]);
    
   for j=1:size(edges)
        x = edges(j, 1);
        y = edges(j, 2);

        for k=1:size(edges)
            % Top left
            if (j ~= k && x == edges(k, 1) - 1 && y == edges(k, 2) - 1)
                hist(1) = hist(1) + 1;
            end

            % Top center
            if (j ~= k && x == edges(k, 1) && y == edges(k, 2) - 1)
                hist(2) = hist(2) + 1;
            end

            % Top right
            if (j ~= k && x == edges(k, 1) + 1 && y == edges(k, 2) - 1)
                hist(3) = hist(3) + 1;
            end

            % Right center
            if (j ~= k && x == edges(k, 1) + 1 && y == edges(k, 2))
                hist(4) = hist(4) + 1;
            end

            % Bottom right
            if (j ~= k && x == edges(k, 1) + 1 && y == edges(k, 2) + 1)
                hist(5) = hist(5) + 1;
            end

            % Bottom center
            if (j ~= k && x == edges(k, 1) && y == edges(k, 2) + 1)
                hist(6) = hist(6) + 1;
            end

            % Bottom left
            if (j ~= k && x == edges(k, 1) - 1 && y == edges(k, 2) + 1)
                hist(7) = hist(7) + 1;
            end

            % Left center
            if (j ~= k && x == edges(k, 1) - 1 && y == edges(k, 2))
                hist(8) = hist(8) + 1;
            end

        end
   end
end