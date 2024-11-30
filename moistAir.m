classdef moistAir < matlab.mixin.Copyable
    %MOISTAIR properties of moist air
    %   _w: water
    %   _da: dry air
    properties(Access=private)
        % These private properties are to avoid the warning of accessing
        % another variable in its setting method. Else without this
        % intermediate level, checking the length of t, rh and q in their
        % setting method cause an warning and unpredicted behaviour
        Priv_t
        Priv_W
        Priv_qm = 1
    end
    properties(Constant,Hidden)    % physical constants
        atm = 101325        % atmospheric pressure in Pa
        MWw = 18.015268     % molecular mass of water in g/mol
        MWda= 28.966        % molecular mass of dry air in g/mol. MWw/MWda = 0.621945
        R   = 8314.472      % universal gas constant
    end
    properties(Dependent)
        t       % temperature in degC
        W       % humidity ratio = mass_w / mass_da in g/g_da
        qm      % mass flow rate kg/s

        Rda     % gas constant of dry air 287.042 J/kgda/K
        Rw      % gas constant of water vapor 461.524 J/kgw/K

        T       % temperature in degK
        rho     % density kg/m3
        rho_da  % density kg_da/m3
        h       % enthalpy kJ/kg
        h_da    % enthalpy kJ/kg_da
        qv      % volume flow rate in m3/s

        pws     % saturate water vapor pressure in Pa at temp T
        xws     % mole frac of water vapor in sat moist air at temp T
        Ws      % humidity ratio = mass_w / mass_da in sat moist air at temp T in g/g_da

        pw      % actual water vapor partial pressure
        xw      % actual mole frac of water vapor at temp T
        rh      % relative humidity 0~1
    end

    methods(Static)
        function pws = p_vsat(t)
            % pws in Pa, T in Kelvin
            T = t + 273.15;
            if t < -100
                error('out of polynomial valid bound');
            elseif t < 0
                C1 = -5.6745359E3;
                C2 =  6.3925247E0;
                C3 = -9.6778430E-3;
                C4 =  6.2215701E-7;
                C5 =  2.0747825E-9;
                C6 = -9.4840240E-13;
                C7 =  4.1635019E0;
                pws = exp(C1./T + C2 + C3*T + C4*T.^2 + C5*T.^3 + C6*T.^4 + C7*log(T));
            elseif t <= 200
                C8 = -5.8002206E3;
                C9 =  1.3914993E0;
                C10= -4.8640239E-2;
                C11=  4.1764768E-5;
                C12= -1.4452093E-8;
                C13=  6.5459673E0;
                pws = exp(C8./T + C9 + C10*T + C11*T.^2 + C12*T.^3 + C13*log(T));
            else
                error('out of polynomial valid bound');
            end
        end
        function [air_in_in, air_out_in] = get_eq_air(air_in_in, air_out_in, ratio)
            % get equilibrium state of moist air in circulation
            % NOTICE: qm is not modified in this function
            air_in = copy(air_in_in); air_out = copy(air_out_in);
            air_mixed = air_in + air_out;
            wind_qm = 1; % kg/s, nominal value
            for i = 1:10
                air_in.W = air_mixed.W;
                air_in.qm  = wind_qm .* (1-ratio);
                air_out.qm = wind_qm .* ratio;
                air_mixed  = air_in + air_out;
                air_mixed  = air_mixed.add_moisture(70/1000/3600,37);
            end

            air_in_in.t = air_in.t;
            air_in_in.W = air_in.W;
            air_out_in.t= air_out.t;
            air_out_in.W= air_out.W;

        end
        function flag = check_accordance(tlen, rhlen, qmlen)
            % to check if length are in accordance
            if (tlen == rhlen) || (tlen == 1) || (rhlen == 1)
                if (tlen == qmlen) || (qmlen == 1)
                    flag = true;
                else
                    error('different length of qm and rh, t!')
                end
            else
                error('different length of t and rh!')
            end
        end
    end

    methods
        function obj = moistAir(tin, rhin)
            %moistAir constructor
            %   check length and set properties
            if moistAir.check_accordance(length(tin), length(rhin), length(obj.qm))
                obj.Priv_t = tin;
                obj.rh     = rhin;
            end
        end

        % set methods
        function set.t(obj,tin)
            moistAir.check_accordance(length(tin), length(obj.W), length(obj.qm));
            if tin > -273.15
                obj.Priv_t = tin;
            else
                error('t in Celsius is less than -273.15! Check input');
            end
        end
        function set.T(obj,Tin)
            % function to synchronize T and t value
            if Tin > 0
                obj.t = Tin - 273.15;
            else
                error('T in Kelvin is less than 0! Check input');
            end
        end
        function set.W(obj, Win)
            moistAir.check_accordance(length(obj.t), length(Win), length(obj.qm));
            obj.Priv_W = Win;
        end
        function set.rh(obj,rhin)
            if rhin<0 | rhin>1
                error('rh is in the range of 0 to 1! Check input')
            else
                pw_t = rhin .* obj.pws;
                obj.Priv_W = obj.MWw / obj.MWda .* pw_t ./ (obj.atm - pw_t);
            end
        end
        function set.h(obj,hin)
            hdain = hin .* (1 + obj.W);
            obj.t = (hdain - 2501 * obj.W) ./ (1.006 + 1.86 * obj.W);
        end
        function set.h_da(obj,hdain)
            hin = hdain./ (1 + obj.W);
            obj.t = (hdain - 2501 * obj.W) ./ (1.006 + 1.86 * obj.W);
        end
        function set.qm(obj,qmin)
            moistAir.check_accordance(length(obj.t), length(obj.W), length(qmin));
            obj.Priv_qm = max(qmin,1e-7);
        end
        function set.qv(obj,qvin)
            obj.qm = qvin .* obj.rho;
        end
        % get methods to update dependent properties
        function t    = get.t(obj)
            t = obj.Priv_t;
        end
        function W    = get.W(obj)
            W = obj.Priv_W;
        end
        function qm    = get.qm(obj)
            qm = obj.Priv_qm;
        end
        function T    = get.T(obj)
            T = obj.t + 273.15;
        end
        function rho  = get.rho(obj)
            rho = obj.rho_da .* (1 + obj.W);
        end
        function rho_da  = get.rho_da(obj)
            % mass_da / volume
            v = obj.Rda * obj.T .* (1 + (1/(obj.MWw/obj.MWda)).*obj.W) / obj.atm;%in ASHRAE Handbook 1.10 (26)
            rho_da = 1./v;
        end
        function h    = get.h(obj)
            h = obj.h_da .* obj.rho_da ./ obj.rho;
        end
        function h_da = get.h_da(obj)
            hda = 1.006 * obj.t;% P1.10 (27)
            hg  = 2501 + 1.86 * obj.t;
            h_da   = hda  + obj.W .* hg;
        end
        function qv   = get.qv(obj)
            qv = obj.qm ./ obj.rho;
        end
        function Rda  = get.Rda(obj)
            Rda = obj.R / obj.MWda;
        end
        function Rw   = get.Rw(obj)
            Rw = obj.R / obj.MWw;
        end
        function pws  = get.pws(obj)
            %pws in Pa, T in Kelvin
            pws = obj.p_vsat(obj.t);
        end
        function xws  = get.xws(obj)
            xws = obj.pws / obj.atm;
        end
        function Ws   = get.Ws(obj)
            Ws = obj.MWw / obj.MWda * obj.pws ./ (obj.atm - obj.pws);
        end
        function pw   = get.pw(obj)
            pw = obj.W * obj.atm ./ (obj.MWw / obj.MWda + obj.W);
        end
        function xw   = get.xw(obj)
            xw = obj.pw / obj.atm;
        end
        function rh   = get.rh(obj)
            rh = obj.pw ./ obj.pws;
        end
        % check saturation and update property
        function obj = update_to_sat(obj)
            % If saturated, the mixture need special treatment on dew/frost
            % Constant latent heat and heat capacity is used.
            % h = h_da + Ws * h_vap + (W-Ws)*h_water,ice
            % h_vap = 2501 + 1.86*t
            % h_water = 4.186*t
            % h_ice = -333.4 + 2.1*t
            % WARNING: No limit is put on input range, be careful with
            % the input temperature!
            if obj.rh > 1
                warning("condense/sublimation occurs, careful with input");
                tol_rel = 1e-3;
                max_iter = 100;
                old_t = 1e6;
                new_h_da = obj.h_da;
                for i = 1:max_iter
                    if abs((obj.t - old_t)/obj.t) < tol_rel
                        break;
                    else
                        if i == max_iter
                            warning("max iter reached, check output for precision")
                        end
                        old_t = obj.t;
                        if obj.t > 0
                            obj.t = (new_h_da -2501 * obj.Ws)...
                                /(1.006 + 1.86 * obj.Ws + 4.186 * (obj.W-obj.Ws));
                        else
                            obj.t = (new_h_da - 2501 * obj.Ws  + 333.4 * obj.W)...
                                /(1.006 + 1.86 * obj.Ws + 2.1 * (obj.W-obj.Ws));
                        end
                        obj.t = (obj.t - old_t)*0.1 + old_t;% relaxation
                    end
                end
                obj.qm = obj.qm.*(1 - (obj.W - obj.Ws)./(1 + obj.W));
                obj.W = obj.Ws;
            end
        end
        % adding moisture to ma (= const P evaporation)
        function obj = add_moisture(obj, qm_w, t_w, varargin)
            % qm_w: mass of water added in kg/s
            % t_w: temperature of water
            % phase: vapor by default
            if nargin == 3
                phase = "v";
            elseif nargin == 4
                phase = varargin{1};
            else
                error("wrong inputs")
            end
            new_W = qm_w ./ (obj.qm ./ (1 + obj.W)) + obj.W;
            if strcmp(phase,"v")
                h_w = 2501 + 1.86 * t_w;
            elseif strcmp(phase,"l")
                h_w = 4.186 * t_w;
            elseif strcmp(phase,"s")
                h_w = -333.4 + 2.1 * t_w;
            end
            new_h_da = (h_w .* qm_w + obj.h .* obj.qm) ./ (qm_w + obj.qm) .* (1 + new_W);
            obj.t = (new_h_da - 2501 * new_W) ./ (1.006 + 1.86 * new_W);
            obj.qm = obj.qm + qm_w;
            obj = obj.update_to_sat();
        end
        % overload "+" to merge two streams
        function obj = plus(obj1, obj2)
            new_qm   = obj1.qm + obj2.qm;
            new_h_da = (obj1.h .* obj1.qm + obj2.h .* obj2.qm) ./ (obj1.qm ./ (1 + obj1.W) + obj2.qm ./ (1 + obj2.W));
            new_W    = (obj1.qm .* (1-(1 ./ (1 + obj1.W))) + obj2.qm .* (1-(1 ./ (1 + obj2.W)))) ...
                ./ (obj1.qm ./ (1 + obj1.W) + obj2.qm ./ (1 + obj2.W));
            new_t    = (new_h_da - 2501*new_W) ./ (1.006 + 1.86 * new_W);
            obj = moistAir(new_t, 0);
            obj.W = new_W;
            obj.qm= new_qm;
            if obj.rh > 1
                obj = obj.update_to_sat();
            end
        end
        % overload "*" to manipulate qm
        function obj = mtimes(obj1, obj2)
            % obj1 is const
            if isa(obj1,'double') && isa(obj2,'moistAir')
                obj = obj2.copy;
                obj.qm = obj.qm * obj1;
                % obj2 is const
            elseif isa(obj1,'moistAir') && isa(obj2,'double')
                obj = obj1.copy;
                obj.qm = obj.qm * obj2;
            end
        end
        % overload subscript indexing
        function sref = subsref(obj, s)
            type = s.type;
            switch type
                case '()'
                    index  = s.subs{1};
                    len_t  = length(obj.t);
                    len_rh = length(obj.rh);
                    len_qm = length(obj.qm);
                    len_max=max([len_t,len_rh,len_qm]);
                    if len_t < len_max
                        tmp_t  = obj.t .* ones(len_max,1);
                    elseif len_rh < len_max
                        tmp_rh = obj.rh .* ones(len_max,1);
                    elseif len_qm < len_max
                        tmp_qm = obj.qm .* ones(len_max,1);
                    end
                    sref = moistAir(obj.t(index), obj.rh(index));
                    sref.qm = obj.qm(index);
                case '{}'
                    sref = builtin('subsref',obj,s);
                case '.'
                    sref = builtin('subsref',obj,s);
            end
        end
    end
end