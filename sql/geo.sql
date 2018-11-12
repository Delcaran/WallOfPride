USE `ios_db`;

--
-- Function elipsoid_P
--
create function elipsoid_P(phi double, e double)
returns double deterministic
return exp(e * atanh(e * sin(phi)));

--
-- Function elipsoid_cos_chi
--
create function elipsoid_cos_chi(phi double, e double)
returns double deterministic
return (2*cos(phi)) / ( ((1+sin(phi)) / elipsoid_P(phi, e)) + ((1-sin(phi))*elipsoid_P(phi, e)) );

--
-- Function elipsoid_sin_chi
--
create function elipsoid_sin_chi(phi double, e double)
returns double deterministic
return ( ((1+sin(phi)) / elipsoid_P(phi, e)) - ((1-sin(phi))*elipsoid_P(phi, e)) ) / ( ((1+sin(phi)) / elipsoid_P(phi, e)) + ((1-sin(phi))*elipsoid_P(phi, e)) );

--
-- Function utm_u
--
create function utm_u(lambda double, phi double, e double)
returns double deterministic
return atanh(elipsoid_cos_chi(phi, e)*sin(lambda));

--
-- Function utm_v
--
create function utm_v(lambda double, phi double, e double)
returns double deterministic
return atan(elipsoid_sin_chi(phi, e) / (elipsoid_cos_chi(phi, e)*cos(lambda)));

--
-- Function utm_f1
--
delimiter //
create function utm_f1(lambda double, phi double, e double, R4 double, a2 double, a4 double, a6 double, a8 double, a10 double, a12 double)
returns double deterministic
begin
    declare u double;
    declare v double;
    set u = utm_u(lambda, phi, e);
    set v = utm_v(lambda, phi, e);
    return R4*(u+ (a2*sinh(2*u)*cos(2*v)) + (a4*sinh(4*u)*cos(4*v)) + (a6*sinh(6*u)*cos(6*v)) + (a8*sinh(8*u)*cos(8*v)) + (a10*sinh(10*u)*cos(10*v)) + (a12*sinh(12*u)*cos(12*v)) );
end//
delimiter ;

--
-- Function utm_f2
--
delimiter //
create function utm_f2(lambda double, phi double, e double, R4 double, a2 double, a4 double, a6 double, a8 double, a10 double, a12 double)
returns double deterministic
begin
    declare u double;
    declare v double;
    set u = utm_u(lambda, phi, e);
    set v = utm_v(lambda, phi, e);
    return R4*(v+ (a2*cosh(2*u)*sin(2*v)) + (a4*cosh(4*u)*sin(4*v)) + (a6*cosh(6*u)*sin(6*v)) + (a8*cosh(8*u)*sin(8*v)) + (a10*cosh(10*u)*sin(10*v)) + (a12*cosh(12*u)*sin(12*v)) );
end//
delimiter ;

--
-- Function utm_x
--
create function utm_x(lambda double, phi double, e double, R4 double, a2 double, a4 double, a6 double, a8 double, a10 double, a12 double, lambda0 double, k0 double, shift double)
returns double deterministic
return k0 * utm_f1(lambda-lambda0, phi, e, R4, a2, a4, a6, a8, a10, a12) + shift;

--
-- Function utm_y
--
create function utm_y(lambda double, phi double, e double, R4 double, a2 double, a4 double, a6 double, a8 double, a10 double, a12 double, lambda0 double, k0 double, shift double)
returns double deterministic
return k0 * utm_f2(lambda-lambda0, phi, e, R4, a2, a4, a6, a8, a10, a12) + shift;

--
-- Function polar_f1
--
create function polar_f1(lambda double, phi double, a double, e double, k90 double)
returns double deterministic
return (2*a*sin(lambda)*elipsoid_cos_chi(phi, e)) / (k90*(1+elipsoid_sin_chi(phi, e)));

--
-- Function polar_f2
--
create function polar_f2(lambda double, phi double, a double, e double, k90 double)
returns double deterministic
return (-2*a*cos(lambda)*elipsoid_cos_chi(phi, e)) / (k90*(1+elipsoid_sin_chi(phi, e)));

--
-- Function ups_x
--
delimiter //
create function ups_x(lambda double, phi double, Z double, k0 double, a double, e double, k90 double, lambda0 double, pole double)
returns double deterministic
begin
    if Z = 1 then
        return k0 * polar_f1(lambda-lambda0, phi, a, e, k90) + pole;
    else
        return k0 * polar_f1(lambda-lambda0, -phi, a, e, k90) + pole;    
    end if;
end//
delimiter ;

--
-- Function ups_y
--
delimiter //
create function ups_y(lambda double, phi double, Z double, k0 double, a double, e double, k90 double, lambda0 double, pole double)
returns double deterministic
begin
    if Z = 1 then
        return k0 * polar_f2(lambda-lambda0, phi, a, e, k90) + pole;
    else
        return -k0 * polar_f2(lambda-lambda0, -phi, a, e, k90) + pole;    
    end if;
end//
delimiter ;


--
-- Function coord_to_DMS
--
delimiter //
CREATE FUNCTION `coord_to_DMS`(coord double, is_lon bool) RETURNS varchar(20) CHARSET cp1251
    NO SQL
    DETERMINISTIC
begin    
  declare acoord double;
  declare sec_coord double;
  declare scoord varchar(20);
  declare padding double;
  
  set padding = 11;
  
  if is_lon then
    set padding = padding + 1;
    if coord>0 then
      set scoord = 'E';
    elseif coord<0 then  
      set scoord = 'W';
    end if;
  else
    if coord>0 then
      set scoord = 'N';
    elseif coord<0 then  
      set scoord = 'S';
    end if;
  end if;
  
  set acoord = abs(coord);
  set sec_coord = acoord*3600;
  set scoord = concat(scoord,' ', lpad(TIME_FORMAT(SEC_TO_TIME(sec_coord), '%H° %i'' %s.'), padding, '0'), rpad(substring_index(cast(sec_coord as char),'.',-1), 3, '0'), '"');
  
  return scoord;
end//
delimiter ;

--
-- Function coord_to_DMM
--
delimiter //
CREATE FUNCTION `coord_to_DMM`(coord double, is_lon bool) RETURNS varchar(20) CHARSET cp1251
    NO SQL
    DETERMINISTIC
begin    
  declare acoord double;
  declare sec_coord double;
  declare scoord varchar(20);
  declare padding double;
  
  set padding = 7;
  
  if is_lon then
    set padding = padding + 1;
    if coord>0 then
      set scoord = 'E';
    elseif coord<0 then  
      set scoord = 'W';
    end if;
  else
    if coord>0 then
      set scoord = 'N';
    elseif coord<0 then  
      set scoord = 'S';
    end if;
  end if;
  
  set acoord = abs(coord);
  set sec_coord = acoord*3600;
  set scoord = concat(scoord,' ', lpad(TIME_FORMAT(SEC_TO_TIME(sec_coord), '%H° %i.'), padding, '0'), rpad(substring_index(cast(sec_coord/60 as char),'.',-1), 3, '0'), '''');
  
  return scoord;
end//
delimiter ;

--
-- Function coord_to_DD
--
delimiter //
CREATE FUNCTION `coord_to_DD`(coord double, is_lon bool) RETURNS varchar(20) CHARSET cp1251
    NO SQL
    DETERMINISTIC
begin    
  declare acoord double;
  declare sec_coord double;
  declare scoord varchar(20);
  declare padding double;
  
  set padding = 2;
  
  if is_lon then
    set padding = padding + 1;
    if coord>0 then
      set scoord = 'E';
    elseif coord<0 then  
      set scoord = 'W';
    end if;
  else
    if coord>0 then
      set scoord = 'N';
    elseif coord<0 then  
      set scoord = 'S';
    end if;
  end if;
  
  set acoord = abs(coord);
  set scoord = concat(scoord,' ', lpad(substring_index(cast(acoord as char),'.',1), padding, '0'), '.', rpad(substring_index(cast(acoord as char),'.',-1), 3, '0'), '°');
  
  return scoord;
end//
delimiter ;

--
-- Function latlon_to_polarUPS
--
delimiter //
CREATE FUNCTION latlon_to_polarUPS(lat double, lon double) RETURNS varchar(200) CHARSET cp1251
    NO SQL
    DETERMINISTIC
begin
    # Constants
    declare a double; # elipsoid semi major axis
    declare b double; # elipsoid semi minor axis
    declare f double; # elipsoid flattening
    declare e double; # elipsoid eccentricity
    declare k90 double;
    declare k0 double; # central scale / point-scale at the Pole, scale at the Pole
    declare lambda0 double; # central meridian / longitude [down|up] from the Pole if [Z=1|Z=-1] (radians)
    declare x_pole double; # easting of the Pole (meters)
    declare y_pole double; # northing of the Pole (meters)    
    # derivates
    declare lambda double; # longitude (radians)
    declare phi double; # latitude (radians)
    declare Z double; # Z=1 if northern emisphere, Z=-1 if southern emisphere
    # outputs
    declare x double; # easting (meters)
    declare y double; # northing (meters)
    
    # CONSTANTS
    
    # WGS84 Elipsoid (in meters)
    set a = 6378137.0;
    set b = 6356752.3142451794976;
    set f = 1 - (b/a);
    set e = sqrt(1 - pow(b/a, 2));
    set k90 = sqrt(1 - pow(e,2)) * exp(e * atanh(e));   
    # Universal Polar Stereographics
    set lambda0 = 0;
    set k0 = 0.994;
    set x_pole = 2000000;
    set y_pole = 2000000;
    
    set lambda = radians(lon);
    set phi = radians(lat);
    
    if lat >= 0 then
        set Z = 1;
    else
        set Z = -1;
    end if;

    set x = ups_x(lambda, phi, Z, k0, a, e, k90, lambda0, x_pole);
    set y = ups_y(lambda, phi, Z, k0, a, e, k90, lambda0, x_pole);
    return concat(cast(x as char),',', cast(y as char));
end//
delimiter ;

--
-- Function latlon_to_UTM
--
delimiter //
CREATE FUNCTION latlon_to_UTM(lat double, lon double) RETURNS varchar(200) CHARSET cp1251
    NO SQL
    DETERMINISTIC
begin
    # Constants
    declare a double; # elipsoid semi major axis
    declare b double; # elipsoid semi minor axis
    declare f double; # elipsoid flattening
    declare e double; # elipsoid eccentricity
    declare k0 double; # central scale
    declare x_cm double; # central meridian easting
    declare y_eq double; # equator northing
    declare R4 double; # meridional isoperimetric radius
    declare a2 double;
    declare a4 double;
    declare a6 double;
    declare a8 double;
    declare a10 double;
    declare a12 double;
    declare phi_origin double;
    declare x_origin double;
    declare y_origin double;
    # derivates
    declare lambda0 double; # central meridian    
    declare lambda_origin double;
    declare lambda double; # longitude (radians)
    declare phi double; # latitude (radians)
    declare Z double; # UTM Zone
    # outputs
    declare x double; # easting (meters)
    declare y double; # northing (meters)
    declare band char(1); # latitude band
	
	if (lat < -80 or lat > 84) then
		return null;
	end if;
    
    # CONSTANTS
    
    # WGS84 Elipsoid (in meters)
    set a = 6378137.0;
    set b = 6356752.3142451794976;
    set f = 1 - (b/a);
    set e = sqrt(1 - pow(b/a, 2));
    set R4 = 6367449.1458234153093;
    set a2 = 8.3773182062446983032 * pow(10, -4);
    set a4 = 7.608527773572489156 * pow(10, -7);
    set a6 = 1.19764550324249210 * pow(10, -9);
    set a8 = 2.4291706803973131 * pow(10, -12);
    set a10 = 5.711818369154105 * pow(10, -15);
    set a12 = 1.47999802705262 * pow(10, -17);
    # UTM
    set lambda = radians(lon);
    set phi = radians(lat);
    set Z = floor((lon+180)/6)+1;
    if phi < 0 then
        set Z = -Z;
    end if;
    if Z = 31 and lat >= 56 and lat < 64 and lon >= 3 then
        set Z = 32;
    elseif Z = 32 and lat >= 72 then
        if lon < 9 then
            set Z = 31;
        else
            set Z = 33;
        end if;
    elseif Z = 34 and lat >= 72 then
        if lon < 21 then
            set Z = 33;
        else
            set Z = 35;
        end if;   
    elseif Z = 36 and lat >= 72 then
        if lon < 33 then
            set Z = 35;
        else
            set Z = 37;
        end if;
    end if;
    
    set lambda0 = radians(-183 + (6 * abs(Z)));
    set lambda_origin = lambda0;
    set phi_origin = 0;
    set k0 = 0.9996;
    set x_origin = 500000;
    if Z>0 then
        set y_origin = 0;
    else
        set y_origin = 10000000;
    end if;
    
    set x_cm = x_origin - k0 * utm_f1(lambda_origin-lambda0, phi_origin, e, R4, a2, a4, a6, a8, a10, a12);
    set y_eq = y_origin - k0 * utm_f2(lambda_origin-lambda0, phi_origin, e, R4, a2, a4, a6, a8, a10, a12);

    set band = substring('CDEFGHJKLMNPQRSTUVWXX', floor(lat/8+10)+1, 1);
    set x = utm_x(lambda, phi, e, R4, a2, a4, a6, a8, a10, a12, lambda0, k0, x_cm);
    set y = utm_y(lambda, phi, e, R4, a2, a4, a6, a8, a10, a12, lambda0, k0, y_eq);    

    return concat(cast(Z as char),',',band,',',cast(x as char),',', cast(y as char));
end//
delimiter ;

--
-- Function print_latlon_to_UTM
--
delimiter //
CREATE FUNCTION print_latlon_to_UTM(lat double, lon double) RETURNS varchar(200) CHARSET cp1251
    NO SQL
    DETERMINISTIC
begin
    declare utm varchar(200);
    declare zone char(2);
    declare band char(1);
    declare x decimal(20,10);
    declare y decimal(20,10);
    set utm = latlon_to_UTM(lat, lon);
    set zone = substring_index(substring_index(utm, ',', 1), ',', -1);
    set band = substring_index(substring_index(utm, ',', 2), ',', -1);
    set x = cast(substring_index(substring_index(utm, ',', 3), ',', -1) as decimal(20,10));
    set y = cast(substring_index(substring_index(utm, ',', 4), ',', -1) as decimal(20,10));
    return concat(zone,' ',band,' ',replace(format(x,0),',',''),' ',replace(format(y,0),',',''));
end//
delimiter ;

--
-- Function latlon_to_MGRS
--
delimiter //
CREATE FUNCTION latlon_to_MGRS(lat double, lon double) RETURNS varchar(200) CHARSET cp1251
    NO SQL
    DETERMINISTIC
begin
    declare utm varchar(200);
    declare ups varchar(200);
    declare zone decimal;
    declare band varchar(1);
    declare x decimal(20,10);
    declare y decimal(20,10);
    declare e100k char(1);
    declare n100k char(1);
    declare z double;
    declare e_ups char(2);
    declare n_ups char(1);
    
    if (lat >= -80 and lat <= 84) then
        set utm = latlon_to_UTM(lat, lon);
        set zone = cast(substring_index(substring_index(utm, ',', 1), ',', -1) as decimal);
        set band = substring_index(substring_index(utm, ',', 2), ',', -1);
        set x = cast(substring_index(substring_index(utm, ',', 3), ',', -1) as decimal(20,10));
        set y = cast(substring_index(substring_index(utm, ',', 4), ',', -1) as decimal(20,10));
        
        if (((zone - 1) % 3) = 0) then 
            set e100k = substring('ABCDEFGH', 1+(floor(x / 100000)-1), 1);
        elseif (((zone - 1) % 3) = 1) then 
            set e100k = substring('JKLMNPQR', 1+(floor(x / 100000)-1), 1);
        else
            set e100k = substring('STUVWXYZ', 1+(floor(x / 100000)-1), 1);
        end if;
        
        if (((zone-1) % 2) = 0) then 
            set n100k = substring('ABCDEFGHJKLMNPQRSTUV', 1+(floor(y / 100000) % 20), 1);
        else
            set n100k = substring('FGHJKLMNPQRSTUVABCDE', 1+(floor(y / 100000) % 20), 1);
        end if;
        
        # tronco a 100km e arrotondo al miglio nautico
        set x = round(x % 100000, 6);
        set y = round(y % 100000, 6);
        
        return concat(cast(zone as char), band, ' ',e100k, n100k, ' ',lpad(replace(format(x,0),',',''), 5,'0'), ' ',lpad(replace(format(y,0),',',''), 5,'0'));
    else
        set ups = latlon_to_polarUPS(lat, lon);
        set x = cast(substring_index(ups, ',', 1) as decimal(20,10));
        set y = cast(substring_index(ups, ',', -1) as decimal(20,10));
        
        if lat >= 0 then
            if (x < 1300000 or x >= 2700000 or y < 1300000 or y >= 2700000) then
                return null;
            end if;
            set e_ups = elt(1+floor(x/100000)-13,'YR','YS','YT','YU','YX','YY','YZ','ZA','ZB','ZC','ZF','ZG','ZH','ZJ');
            set n_ups = elt(1+floor(y/100000)-13,'A','B','C','D','E','F','G','H','J','K','L','M','N','P');
        else
            if (x < 800000 or x >= 3200000 or y < 800000 or y >= 3200000) then
                return null;
            end if;
            set e_ups = elt(1+floor(x/100000)-8,'AJ','AK','AL','AP','AQ','AR','AS','AT','AU','AX','AY','AZ','BA','BB','BC','BF','BG','BH','BJ','BK','BL','BP','BQ','BR');
            set n_ups = elt(1+floor(y/100000)-8,'A','B','C','D','E','F','G','H','J','K','L','M','N','P','Q','R','S','T','U','V','W','X','Y','Z');
        end if;
        set x = round(x % 100000, 6);
        set y = round(y % 100000, 6);
        return concat(e_ups, n_ups, ' ',lpad(replace(format(x,0),',',''), 5,'0'), ' ',lpad(replace(format(y,0),',',''), 5,'0'));
    end if;
end//
delimiter ;
