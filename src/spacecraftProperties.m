function sc = spacecraftProperties(bus,solar,onePanelFlag)
%SPACECRAFTPROPERTIES Computes spacecraft properties (CG, MoI, etc)
%   Given a bus parameters (dimensions, cg, MoI, etc) and solar panel
%   properties, computes the resulting spacecraft's properties

% setup
arguments
    bus
    solar
    onePanelFlag = false
end

M2B_rot = rotz(90)';

% mass
sc.mass = bus.mass + 2*solar.mass;

% center of gravity/mass (cg)
rBus2SolarLeft_M = [-solar.offsetInX_M - solar.dim.x_M, solar.axisXLocation_M(1), solar.axisXLocation_M(2)];
rBus2SolarRight_M = [bus.dim.x_M + solar.offsetInX_M, solar.axisXLocation_M(1), solar.axisXLocation_M(2)];

sc.cg.x_M = (bus.mass*bus.cg.x_M + solar.mass*(rBus2SolarLeft_M(1) + solar.cg.x_M) + solar.mass*(rBus2SolarRight_M(1) + solar.cg.x_M))/(bus.mass + 2*solar.mass);
sc.cg.y_M = (bus.mass*bus.cg.y_M + solar.mass*(rBus2SolarLeft_M(2) + solar.cg.y_M) + solar.mass*(rBus2SolarRight_M(2) + solar.cg.y_M))/(bus.mass + 2*solar.mass);
sc.cg.z_M = (bus.mass*bus.cg.z_M + solar.mass*(rBus2SolarLeft_M(3) + solar.cg.z_M) + solar.mass*(rBus2SolarRight_M(3) + solar.cg.z_M))/(bus.mass + 2*solar.mass);

% moment of inertia (MoI)
cgOffsetBus_B = M2B_rot*([bus.cg.x_M, bus.cg.y_M, bus.cg.z_M] - [sc.cg.x_M, sc.cg.y_M, sc.cg.z_M])';
cgOffsetSolarLeft_B = M2B_rot*(rBus2SolarLeft_M + [solar.cg.x_M, solar.cg.y_M, solar.cg.z_M] - [sc.cg.x_M, sc.cg.y_M, sc.cg.z_M])';
cgOffsetSolarRight_B = M2B_rot*(rBus2SolarRight_M + [solar.cg.x_M, solar.cg.y_M, solar.cg.z_M] - [sc.cg.x_M, sc.cg.y_M, sc.cg.z_M])';

sc.MoI.x_B = bus.MoI.x_B + bus.mass*norm(cgOffsetBus_B(2:3))^2 + 2*solar.MoI.x_B + ...
    solar.mass*norm(cgOffsetSolarLeft_B(2:3))^2 + solar.mass*norm(cgOffsetSolarRight_B(2:3))^2;
sc.MoI.y_B = bus.MoI.y_B + bus.mass*norm(cgOffsetBus_B([1 3]))^2 + 2*solar.MoI.x_B + ...
    solar.mass*norm(cgOffsetSolarLeft_B([1 3]))^2 + solar.mass*norm(cgOffsetSolarRight_B([1 3]))^2;
sc.MoI.z_B = bus.MoI.z_B + bus.mass*norm(cgOffsetBus_B(1:2))^2 + 2*solar.MoI.x_B + ...
    solar.mass*norm(cgOffsetSolarLeft_B(1:2))^2 + solar.mass*norm(cgOffsetSolarRight_B(1:2))^2;

end

