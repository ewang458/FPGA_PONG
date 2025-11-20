library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity debounce is
	port(
		clk:   in  std_logic;
		btn: in  std_logic;
		btn_out: out std_logic
	);
end debounce;

architecture arch of debounce is


	signal cnt : unsigned (21 downto 0):= to_unsigned(2500000,22);
	signal btn_sr: std_logic_vector(2 downto 0);
	signal btnState: std_logic := '0';
begin

btn_out <= btnState;
process(clk)
begin
if rising_edge(clk) then
    btn_sr <= btn_sr(1 downto 0) & btn;

-- Saturation counter
    if (btn_sr(2)= '0') then
        if (cnt/= 0) then
            cnt<=cnt-1;
        else
            btnState <= '1';
        end if;
    elsif (btn_sr(2)= '1') then
        if (cnt<2500000) then
            cnt <= cnt + 1;
        else 
            btnState <= '0';
        end if;
    end if;
    
end if;
end process;
end arch;