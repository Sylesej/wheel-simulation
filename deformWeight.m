function newSpokes = deformWeight(spokes,fWeight,wheeldata)
%This function calculates deformation of a wheel given the force and data
%about the wheel, and outputs changed spoke coordinates.
%Inputs:
%spokes is an array of coordinates of spoke endpoints.
%fWeight is a number specifying the amplitude of force. Assumed
%perpendicular to the rim.
%wheeldata is a struct containing data such as bending stiffness.

newSpokes = spokes;