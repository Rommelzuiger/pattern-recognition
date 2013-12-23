function im = owndata2im(data)
    im = zeros([16 16]);
    for i=1:16
        im(i,:) = data(1 + (i-1)*16:i*16);
    end
    im = im';
end
        