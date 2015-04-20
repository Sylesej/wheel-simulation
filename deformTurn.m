function newSpokes = deformTurn(spokes,fTurn,wheeldata)
%This function calculates deformation of a wheel given the force and data
%about the wheel, and outputs changed spoke coordinates.
%Inputs:
%spokes is an array of coordinates of spoke endpoints.
%fTurn is a number specifying the amplitude of force. Assumed
%perpendicular to the rim plane and driving direction
%wheeldata is a struct containing data such as bending stiffness.

newSpokes = spokes;