function [outputRefVec,outputStateEstVec] = demux(outputRefTrajectory, outputStateEst)
%DEMUX Summary of this function goes here
%   Detailed explanation goes here
el = [1 4; 2 5; 3 6];
for i = 1:3
    outputRefVec(i).xRefI = outputRefTrajectory.xRefI(el(i,:),:);
    outputRefVec(i).vRefI = outputRefTrajectory.vRefI(i,:);
    outputRefVec(i).xRefT = outputRefTrajectory.xRefT(el(i,:),:);
    outputRefVec(i).vRefT = outputRefTrajectory.vRefT(i,:);
    outputStateEstVec(i).chaserState = outputStateEst.chaserState(el(i,:),1);
end
end

