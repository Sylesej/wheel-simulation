function plotWheel(spokes,deform)
%This function plots the wheel, given the spoke coordinates given in the
%spokes array.
nSpokes = length(spokes(:,1));

plotSpokeX = zeros(1,2*nSpokes);
plotSpokeY = zeros(1,2*nSpokes);
plotSpokeZ = zeros(1,2*nSpokes);

for n = 1:nSpokes
    plotSpokeX(n*2-1) = spokes(n,1);
    plotSpokeX(n*2)   = spokes(n,2);
    plotSpokeY(n*2-1) = spokes(n,3);
    plotSpokeY(n*2)   = spokes(n,4);
    plotSpokeZ(n*2-1) = spokes(n,5);
    plotSpokeZ(n*2)   = spokes(n,6); 
end

for n = 1:2:2*nSpokes
    plot3([plotSpokeX(n) plotSpokeX(n+1)],[plotSpokeY(n) plotSpokeY(n+1)],...
        [plotSpokeZ(n) plotSpokeZ(n+1)],'r')
    if n==1
        hold on
    end
end
%plot3(plotSpoke(n*2-1),plotSpokeY,plotSpokeZ)
axis 'equal'
xlabel('z [m]')
ylabel('y [m]')
zlabel('x [m]')
if deform==0
    title('Udeformeret hjul')
else
    title('Deformeret hjul')
end
plot3(spokes(:,2),spokes(:,4),spokes(:,6),'g')
%legend('Spokes','Rim')