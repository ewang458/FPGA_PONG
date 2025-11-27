library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;
library UNISIM;
use UNISIM.vcomponents.all;

entity score is
    generic(
        MAX_SCORE : natural := 11
    );
   	port(
		clk:   in    std_logic;
        score_flag  : in  std_logic;
        identity: in std_logic;
        hcount:   in unsigned(9 downto 0);
        vcount:   in unsigned(9 downto 0);
        char_pixel : out std_logic;
        current_score : out integer range 0 to 20;
        rst : in std_logic
	);
end score;

architecture arch of score is
    signal digit0, digit1: integer range 0 to 9; -- digit0->ones digit1->tens
    signal score_internal : integer range 0 to 20 := 0;
    signal prev_flag: std_logic := '0';
    
    -- Font ROM (8x8 patterns for digits 0-9)
    type font_rom_type is array (0 to 9, 0 to 7) of std_logic_vector(7 downto 0);
    constant FONT_ROM : font_rom_type := (
        -- 0
        0 => (X"3C", X"66", X"66", X"66", X"66", X"66", X"66", X"3C"),
        -- 1
        1 => (X"18", X"38", X"18", X"18", X"18", X"18", X"18", X"7E"),
        -- 2
        2 => (X"3C", X"66", X"06", X"0C", X"18", X"30", X"60", X"7E"),
        -- 3
        3 => (X"3C", X"66", X"06", X"1C", X"06", X"06", X"66", X"3C"),
        -- 4
        4 => (X"0C", X"1C", X"2C", X"4C", X"7E", X"0C", X"0C", X"0C"),
        -- 5
        5 => (X"7E", X"60", X"60", X"7C", X"06", X"06", X"66", X"3C"),
        -- 6
        6 => (X"3C", X"66", X"60", X"7C", X"66", X"66", X"66", X"3C"),
        -- 7
        7 => (X"7E", X"06", X"0C", X"18", X"30", X"30", X"30", X"30"),
        -- 8
        8 => (X"3C", X"66", X"66", X"3C", X"66", X"66", X"66", X"3C"),
        -- 9
        9 => (X"3C", X"66", X"66", X"66", X"3E", X"06", X"66", X"3C")
    );
    
begin

    process(clk)
    begin
        if rising_edge(clk) then
            --check for game reset
            if rst = '1' then
                score_internal <= 0;
                prev_flag <= '0';
            else
                if score_flag = '1' and prev_flag = '0' then         -- Detect rising edge of score_flag
                    if score_internal < MAX_SCORE then
                        score_internal <= score_internal + 1;
                    end if;
                end if;
                prev_flag <= score_flag;
            end if;
        end if;
    end process;

    current_score <= score_internal;
    --score integer to digits
    process(score_internal)
        variable temp: integer range 0 to 20;
    begin
        temp := score_internal;
        digit0 <= temp mod 10; --ones digit
        digit1 <= (temp / 10) mod 10; --tens digit
    end process;
    
    --render digits on vga:
    process(hcount, vcount, digit0, digit1, identity)
        variable digit_index: integer range 0 to 1;
        variable current_digit: integer range 0 to 9;
        variable char_start_x: integer;
        variable char_start_y: integer;
        variable x_in_char, y_in_char: integer;
        variable font_x, font_y: integer range 0 to 7;
    begin
        char_pixel <= '0';
        
        if identity = '0' then
            char_start_x := 120;
        else
            char_start_x := 480;
        end if;
        char_start_y := 10;
        
        -- For single digit (0-9), only show ones digit
        -- For double digit (10-20), show both digits
        if digit1 = 0 then
            -- Single digit: only display at position 0
            if to_integer(hcount) >= char_start_x and 
               to_integer(hcount) < char_start_x + 24 and
               to_integer(vcount) >= char_start_y and 
               to_integer(vcount) < char_start_y + 24 then
                
                x_in_char := to_integer(hcount) - char_start_x;
                y_in_char := to_integer(vcount) - char_start_y;
                font_x := x_in_char / 3;
                font_y := y_in_char / 3;
                current_digit := digit0;
                
                if FONT_ROM(current_digit, font_y)(7 - font_x) = '1' then
                    char_pixel <= '1';
                end if;
            end if;
        else
            -- Double digit: display both
            if to_integer(hcount) >= char_start_x and 
               to_integer(hcount) < char_start_x + (2 * 24) and
               to_integer(vcount) >= char_start_y and 
               to_integer(vcount) < char_start_y + 24 then
                
                digit_index := (to_integer(hcount) - char_start_x) / 24;
                x_in_char := (to_integer(hcount) - char_start_x) mod 24;
                y_in_char := to_integer(vcount) - char_start_y;
                font_x := x_in_char / 3;
                font_y := y_in_char / 3;
                
                case digit_index is
                    when 0 => current_digit := digit1;
                    when 1 => current_digit := digit0;
                    when others => current_digit := 0;
                end case;
                
                if FONT_ROM(current_digit, font_y)(7 - font_x) = '1' then
                    char_pixel <= '1';
                end if;
            end if;
        end if;
    end process;
    
end arch;
