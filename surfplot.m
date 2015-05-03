function surfplot(data)

surf(data);
daspect([1,1,(max(data(:))-min(data(:)))/100])
axis tight
shading flat

end