library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;
library UNISIM;
use UNISIM.vcomponents.all;

entity lab_f is
	port(
		clk:   in    std_logic;
		tx:    out   std_logic;
		red:   out   std_logic_vector(1 downto 0);
		green: out   std_logic_vector(1 downto 0);
		blue:  out   std_logic_vector(1 downto 0);
		hsync: out   std_logic;
		vsync: out   std_logic;
		up1: in std_logic;
		up2: in std_logic;
		down1: in std_logic;
		down2: in std_logic;
		reset : in std_logic;
		start :in std_logic
	);
end lab_f;

architecture arch of lab_f is

component debounce 
    port(
        clk: in std_logic;
        btn: in std_logic;
        btn_out: out std_logic
        );
end component;   
	signal clkfb:    std_logic;
	signal clkfx:    std_logic;
	signal hcount:   unsigned(9 downto 0);
	signal vcount:   unsigned(9 downto 0);
	signal hcount_d: unsigned(9 downto 0);
	signal vcount_d: unsigned(9 downto 0);
	signal blank:    std_logic;
	signal frame:    std_logic;
	signal grid_val: std_logic_vector(1 downto 0);
	signal pad2_blu: std_logic_vector(1 downto 0);
	signal pad2_top: unsigned(5 downto 0) := to_unsigned(23, 6);
	signal pad2_bot: unsigned(5 downto 0) := to_unsigned(35, 6);
	signal pad1_top: unsigned(5 downto 0) := to_unsigned(23, 6);
	signal pad1_bot: unsigned(5 downto 0) := to_unsigned(35, 6);
	--signal pad2_top: unsigned(4 downto 0) := to_unsigned(11, 5);
	--signal pad2_bot: unsigned(4 downto 0) := to_unsigned(17, 5);
	--signal pad1_top: unsigned(4 downto 0) := to_unsigned(11, 5);
	--signal pad1_bot: unsigned(4 downto 0) := to_unsigned(17, 5);
	signal ball_dx: integer:= 1;
	signal ball_dy: integer:= 1;
	signal y_change: integer; 
	signal y_thresh: integer;
	signal y_count: unsigned(3 downto 0);
	signal ball_row: unsigned(5 downto 0):= to_unsigned(30,6);
	signal ball_col: unsigned(6 downto 0) := to_unsigned(40,7);
	--signal ball_row: unsigned(4 downto 0):= to_unsigned(15,5);
	--signal ball_col: unsigned(5 downto 0) := to_unsigned(20,6);
	signal ballcnt: unsigned(2 downto 0);
	signal score1: std_logic:= '0';
	signal score2: std_logic:= '0';
	signal hit_counter: unsigned(2 downto 0);
	signal speed: unsigned(2 downto 0) := to_unsigned(4,3);
	signal up_1: std_logic;
	signal up_2: std_logic;
	signal down_1: std_logic;
	signal down_2: std_logic;
	signal rst: std_logic;
	signal h_count_d1: unsigned(9 downto 0);
	signal h_count_d2: unsigned(9 downto 0);
	signal h_count_d3: unsigned(9 downto 0);
	type bg is array(59 downto 0,79  downto 0) of std_logic_vector(1 downto 0);
	--type bg is array(29 downto 0,39  downto 0) of std_logic_vector(1 downto 0);
	signal pong_bg: bg;
	signal xcount: std_logic := '0';
	signal rowcnt: unsigned(5 downto 0) := (others => '0');
	signal colcnt : unsigned(6 downto 0) := (others => '0');
	--signal rowcnt: unsigned(4 downto 0) := (others => '0');
	--signal colcnt : unsigned(5 downto 0) := (others => '0');
	type FSM is (startFrame, startGame, idle, initialize, updatePaddle, updateBall, moveBall, addPaddle, addBall);
	signal game: FSM;

begin
    up1debounce: debounce port map(clk => clkfx, btn => up1, btn_out => up_1);
    up2debounce: debounce port map(clk => clkfx, btn => up2, btn_out => up_2);
    down1debounce: debounce port map(clk => clkfx, btn => down1, btn_out => down_1);
    down2debounce: debounce port map(clk => clkfx, btn => down2, btn_out => down_2);
	tx<='1';

	------------------------------------------------------------------
	-- Clock management tile
	--
	-- Input clock: 12 MHz
	-- Output clock: 25.2 MHz
	--
	-- CLKFBOUT_MULT_F: 50.875
	-- CLKOUT0_DIVIDE_F: 24.250
	-- DIVCLK_DIVIDE: 1
	------------------------------------------------------------------
	cmt: MMCME2_BASE generic map (
		-- Jitter programming (OPTIMIZED, HIGH, LOW)
		BANDWIDTH=>"OPTIMIZED",
		-- Multiply value for all CLKOUT (2.000-64.000).
		CLKFBOUT_MULT_F=>50.875,
		-- Phase offset in degrees of CLKFB (-360.000-360.000).
		CLKFBOUT_PHASE=>0.0,
		-- Input clock period in ns to ps resolution (i.e. 33.333 is 30 MHz).
		CLKIN1_PERIOD=>83.333,
		-- Divide amount for each CLKOUT (1-128)
		CLKOUT1_DIVIDE=>1,
		CLKOUT2_DIVIDE=>1,
		CLKOUT3_DIVIDE=>1,
		CLKOUT4_DIVIDE=>1,
		CLKOUT5_DIVIDE=>1,
		CLKOUT6_DIVIDE=>1,
		-- Divide amount for CLKOUT0 (1.000-128.000):
		CLKOUT0_DIVIDE_F=>24.250,
		-- Duty cycle for each CLKOUT (0.01-0.99):
		CLKOUT0_DUTY_CYCLE=>0.5,
		CLKOUT1_DUTY_CYCLE=>0.5,
		CLKOUT2_DUTY_CYCLE=>0.5,
		CLKOUT3_DUTY_CYCLE=>0.5,
		CLKOUT4_DUTY_CYCLE=>0.5,
		CLKOUT5_DUTY_CYCLE=>0.5,
		CLKOUT6_DUTY_CYCLE=>0.5,
		-- Phase offset for each CLKOUT (-360.000-360.000):
		CLKOUT0_PHASE=>0.0,
		CLKOUT1_PHASE=>0.0,
		CLKOUT2_PHASE=>0.0,
		CLKOUT3_PHASE=>0.0,
		CLKOUT4_PHASE=>0.0,
		CLKOUT5_PHASE=>0.0,
		CLKOUT6_PHASE=>0.0,
		-- Cascade CLKOUT4 counter with CLKOUT6 (FALSE, TRUE)
		CLKOUT4_CASCADE=>FALSE,
		-- Master division value (1-106)
		DIVCLK_DIVIDE=>1,
		-- Reference input jitter in UI (0.000-0.999).
		REF_JITTER1=>0.0,
		-- Delays DONE until MMCM is locked (FALSE, TRUE)
		STARTUP_WAIT=>FALSE
	) port map (
		-- User Configurable Clock Outputs:
		CLKOUT0=>clkfx,  -- 1-bit output: CLKOUT0
		CLKOUT0B=>open,  -- 1-bit output: Inverted CLKOUT0
		CLKOUT1=>open,   -- 1-bit output: CLKOUT1
		CLKOUT1B=>open,  -- 1-bit output: Inverted CLKOUT1
		CLKOUT2=>open,   -- 1-bit output: CLKOUT2
		CLKOUT2B=>open,  -- 1-bit output: Inverted CLKOUT2
		CLKOUT3=>open,   -- 1-bit output: CLKOUT3
		CLKOUT3B=>open,  -- 1-bit output: Inverted CLKOUT3
		CLKOUT4=>open,   -- 1-bit output: CLKOUT4
		CLKOUT5=>open,   -- 1-bit output: CLKOUT5
		CLKOUT6=>open,   -- 1-bit output: CLKOUT6
		-- Clock Feedback Output Ports:
		CLKFBOUT=>clkfb,-- 1-bit output: Feedback clock
		CLKFBOUTB=>open, -- 1-bit output: Inverted CLKFBOUT
		-- MMCM Status Ports:
		LOCKED=>open,    -- 1-bit output: LOCK
		-- Clock Input:
		CLKIN1=>clk,   -- 1-bit input: Clock
		-- MMCM Control Ports:
		PWRDWN=>'0',     -- 1-bit input: Power-down
		RST=>'0',        -- 1-bit input: Reset
		-- Clock Feedback Input Port:
		CLKFBIN=>clkfb  -- 1-bit input: Feedback clock
	);

	------------------------------------------------------------------
	-- VGA display counters
	--
	-- Pixel clock: 25.175 MHz (actual: 25.2 MHz)
	-- Horizontal count (active low sync):
	--     0 to 639: Active video
	--     640 to 799: Horizontal blank
	--     656 to 751: Horizontal sync (active low)
	-- Vertical count (active low sync):
	--     0 to 479: Active video
	--     480 to 524: Vertical blank
	--     490 to 491: Vertical sync (active low)
	------------------------------------------------------------------
	process(clkfx)
	begin
		if rising_edge(clkfx) then
			-- Pixel position counters
			if (hcount>=to_unsigned(799,10)) then
				hcount<=(others=>'0');
				if (vcount>=to_unsigned(524,10)) then
					vcount<=(others=>'0');
				else
					vcount<=vcount+1;
				end if;
			else
				hcount<=hcount+1;
			end if;
			-- Sync, blank and frame
			if (hcount>=to_unsigned(656,10)) and
				(hcount<=to_unsigned(751,10)) then
				hsync<='0';
			else
				hsync<='1';
			end if;
			if (vcount>=to_unsigned(490,10)) and
				(vcount<=to_unsigned(491,10)) then
				vsync<='0';
			else
				vsync<='1';
			end if;
			if (hcount>=to_unsigned(640,10)) or
				(vcount>=to_unsigned(480,10)) then
				blank<='1';
			else
				blank<='0';
			end if;
			if (hcount=to_unsigned(639,10)) and
				(vcount=to_unsigned(479,10)) then
				frame<='1';
			else
				frame<='0';
			end if;
		end if;
	end process;

	------------------------------------------------------------------
	-- VGA output with blanking
	------------------------------------------------------------------
	red<=b"00" when blank='1' else pong_bg(to_integer(vcount(9 downto 3)), to_integer(hcount_d(9 downto 3))); --downto3
	green<=b"00" when blank='1' else pong_bg(to_integer(vcount(9 downto 3)), to_integer(hcount_d(9 downto 3)));
	blue<=b"00" when blank='1' else pong_bg(to_integer(vcount(9 downto 3)), to_integer(hcount_d(9 downto 3)));
	
	process(clkfx) 
	begin 
	   if (rising_edge(clkfx)) then 
	       hcount_d <= hcount;
	   end if;
	end process;
    
    process(clkfx)
    variable ball_y:integer;
    begin
    if (reset = '1') then
        ball_row <= to_unsigned(30,6);
        ball_col <= to_unsigned(40,7);
        pad1_top <= to_unsigned(23,6);
        pad1_bot <= to_unsigned(35,6);
        pad2_top <= to_unsigned(23,6);
        pad2_bot <= to_unsigned(35,6);
        --ball_row <= to_unsigned(15,5);
        --ball_col <= to_unsigned(20,6);
        --pad1_top <= to_unsigned(11,5);
        --pad1_bot <= to_unsigned(17,5);
        --pad2_top <= to_unsigned(11,5);
        --pad2_bot <= to_unsigned(17,5);
        ball_dx <= 1;
        ball_dy <= 1;
        score1 <= '0';
        score2 <= '0';
        xcount <= '0';
        speed <= to_unsigned(4,3);
        ballcnt <= to_unsigned(0,3);
        hit_counter <= to_unsigned(0,3);
        rowcnt <= to_unsigned(0,6);
        colcnt <= to_unsigned(0,7);
        --rowcnt <= to_unsigned(0,5);
        --colcnt <= to_unsigned(0,6);
        game <= startFrame;
        
    elsif (rising_edge(clkfx)) then
        case (game) is --add game start and game end states and reset capabilty 
            when startFrame => 
                ball_row <= to_unsigned(30,6);
                ball_col <= to_unsigned(40, 7);
                if (rowcnt < 60) then 
                    if (colcnt < 80) then 
                        if ((colcnt = 0) or (colcnt = 79)) then 
                            if ((rowcnt >= pad1_top) and (rowcnt < pad1_bot)) then 
                                pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"11";
                            else
                                if ((rowcnt > 0) and (rowcnt < 59)) then
                                    pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"00";
                                elsif ((rowcnt = 0) or (rowcnt = 59)) then 
                                    pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"11";
                                end if;
                            end if;
                        elsif ((colcnt = 39) or (colcnt = 40)) then 
                            pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"11";
                        elsif ((rowcnt = 0) or (rowcnt = 59)) then 
                            pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"11";
                        else                     
                            pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"00";
                        end if;
                        colcnt <= colcnt + 1;
                    else 
                        colcnt <= to_unsigned(0,7);
                        rowcnt <= rowcnt + 1;
                    end if;
                else 
                    rowcnt <= to_unsigned(0,6);
                    game <= startGame;
                end if;                      
            when startGame =>
                hit_counter <= to_unsigned(0,3);
                speed <= to_unsigned(4,3);
                score1 <= '0';
                score2 <= '0';
                if (start = '1') then 
                    game <= idle;
                end if;
            when idle => 
                if (frame = '1') then 
                    game <= initialize;
                end if;
            when initialize => 
                -- score here 
                --60, 39, 40
                if (rowcnt < 60) then
                    pong_bg(to_integer(rowcnt), 39) <= b"11";
                    pong_bg(to_integer(rowcnt), 40) <= b"11";
                    rowcnt <= rowcnt + 1;
                else --80, 0, 59
                    if (colcnt < 80) then 
                        pong_bg(0, to_integer(colcnt)) <= b"11";
                        pong_bg(59, to_integer(colcnt)) <= b"11";
                        colcnt <= colcnt + 1;
                    else --7, 6
                        colcnt <= to_unsigned(0,7);
                        rowcnt <= to_unsigned(0,6);
                        game <= updatePaddle;
                    end if;
                end if;  
            when updatePaddle =>
                if (up_1 = '1')then
                   if (pad1_top > 1) then 
                       pad1_top <= pad1_top - 1;
                       pad1_bot <= pad1_bot - 1;
                   end if;
	            elsif (down_1 = '1') then --59 
	               if (pad1_bot < 59) then
                       pad1_top <= pad1_top + 1;
                       pad1_bot <= pad1_bot + 1;
                   end if;
                end if; 
               if (up_2 = '1')then
                   if (pad2_top > 1) then 
                       pad2_top <= pad2_top - 1;
                       pad2_bot <= pad2_bot - 1;
                   end if;
	            elsif (down_2 = '1') then 
	               if (pad2_bot < 59) then --59
                       pad2_top <= pad2_top + 1;
                       pad2_bot <= pad2_bot + 1;
                   end if;
                end if; 
                game <= updateBall; -- add update score state 
            when updateBall =>      -- if in update score block then only set to 1 if ball, no setting to zero otherwise          
                        if ((ball_col > 1) and (ball_col < 78)) then --1, 78
                            if ((ball_row <= 1) or (ball_row >= 58))then 
                                ball_y := -ball_dy;
                            else 
                                ball_y := ball_dy;
                            end if;
                        elsif ((ball_col = 1) or (ball_col = 78)) then --1, 78
                            if (ball_col = 1) then 
                                if ((ball_row >= pad1_top) and (ball_row < pad1_top + 3)) then  -- make sure u adjust to 12 
                                   ball_col <= to_unsigned(2,7);
                                    if (ball_dy < 0) then 
                                        if ((ball_dy - 1) < -4) then 
                                            ball_y := -4;
                                            ball_dy <= -4;
                                        else
                                            ball_y:= ball_dy -1;
                                            ball_dy <= ball_dy -1;
                                        end if;
                                    else  
                                        if ((-ball_dy - 1) < -4) then
                                            ball_y := -4;
                                            ball_dy <= -4;
                                        else 
                                            ball_y := -ball_dy -1;
                                            ball_dy <= -ball_dy -1;
                                        end if;
                                    end if;
                                    ball_dx <= -ball_dx;
                                    hit_counter <= hit_counter + 1;
                                elsif ((ball_row >= pad1_top +3) and (ball_row <pad1_top +9)) then 
                                  ball_col <= to_unsigned(2,7);
                                    if (ball_dy < 0) then 
                                        ball_y := -ball_dy - 1; 
                                        ball_dy <= -ball_dy -1;  
                                    else 
                                        ball_y := -ball_dy + 1;
                                        ball_dy <= -ball_dy +1;
                                    end if;
                                    ball_dx <= -ball_dx;
                                    hit_counter <= hit_counter + 1;
                                elsif ((ball_row >= pad1_top + 9) and (ball_row < pad1_bot)) then 
                                   ball_col <= to_unsigned(2,7);
                                    if (ball_dy < 0) then 
                                        if (-ball_dy + 1 > 4) then 
                                            ball_y := 4;
                                            ball_dy <= 4;
                                        else 
                                            ball_dy <= -ball_dy + 1;
                                            ball_y := ball_dy + 1;
                                        end if;
                                    elsif (ball_dy > 0) then 
                                        if (ball_dy + 1 > 4) then 
                                            ball_y := 4;
                                            ball_dy <= 4;
                                        else 
                                            ball_dy <= ball_dy + 1;
                                            ball_y := ball_dy + 1;
                                        end if;
                                    else
                                        ball_y := 1;
                                        ball_dy <= 1;
                                    end if;
                                    ball_dx <= -ball_dx;
                                    hit_counter <= hit_counter + 1;
                                else
                                    ball_y:= ball_dy;
                                end if;
                            elsif (ball_col = 78) then --78
                                if ((ball_row >= pad2_top) and (ball_row < pad2_top + 3)) then  -- make sure u adjust to 12 
                                   ball_col <= to_unsigned(77,7);
                                    if (ball_dy < 0) then 
                                        if ((ball_dy - 1) < -4) then 
                                            ball_y := -4;
                                            ball_dy <= -4;
                                        else
                                            ball_y:= ball_dy -1;
                                            ball_dy <= ball_dy -1;
                                        end if;
                                    else  
                                        if ((-ball_dy - 1) < -4) then
                                            ball_y := -4;
                                            ball_dy <= -4;
                                        else 
                                            ball_y := -ball_dy -1;
                                            ball_dy <= -ball_dy -1;
                                        end if;
                                    end if;
                                    ball_dx <= -ball_dx;
                                    hit_counter <= hit_counter + 1;
                                elsif ((ball_row >= pad2_top +3) and (ball_row <pad2_top +9)) then 
                                   ball_col <= to_unsigned(77,7);
                                    if (ball_dy < 0) then 
                                        ball_y := -ball_dy - 1; 
                                        ball_dy <= -ball_dy -1;  
                                    else 
                                        ball_y := -ball_dy + 1;
                                        ball_dy <= -ball_dy +1;
                                    end if;
                                    ball_dx <= -ball_dx;
                                    hit_counter <= hit_counter + 1;
                                elsif ((ball_row >= pad2_top + 9) and (ball_row < pad2_bot)) then 
                                   ball_col <= to_unsigned(77,7);
                                    if (ball_dy < 0) then 
                                        if (-ball_dy + 1 > 4) then 
                                            ball_y := 4;
                                            ball_dy <= 4;
                                        else 
                                            ball_dy <= -ball_dy + 1;
                                            ball_y := ball_dy + 1;
                                        end if;
                                    elsif (ball_dy > 0) then 
                                        if (ball_dy + 1 > 4) then 
                                            ball_y := 4;
                                            ball_dy <= 4;
                                        else 
                                            ball_dy <= ball_dy + 1;
                                            ball_y := ball_dy + 1;
                                        end if;
                                    else
                                        ball_y := 1;
                                        ball_dy<= 1;
                                    end if;
                                    ball_dx <= -ball_dx;
                                    hit_counter <= hit_counter + 1;
                                else
                                    ball_y:= ball_dy;
                                end if;
                          end if;
                      elsif ((ball_col = 0) or (ball_col = 79)) then  -- --0 or 79can edit to go to middle and then go again when start button hit from users side but for now this good
                            if (ball_col = 0) then
                                ball_dx <= 1;
                                if (ball_dy < 0) then 
                                    ball_y:= -1;
                                    ball_dy <= -1;
                                else 
                                    ball_y:= 1;
                                    ball_dy <= 1;
                                end if;
                                score2<= '1';
                             elsif (ball_col = 79) then --79
                                ball_dx <= -1;
                                if (ball_dy < 0) then 
                                    ball_y := -1;
                                    ball_dy <= -1;
                                else 
                                    ball_y := 1;
                                    ball_dy <= 1;
                                end if;
                                score1 <= '1';
                             end if;
                    end if;
                    if (ball_dy /= ball_y) then 
                        y_count <= to_unsigned(0,4);
                        if (ball_y < 0) then 
                            y_thresh <= -ball_y;
                            y_change <= -1;
                        else 
                            y_thresh <= ball_y;
                            y_change <= 1;
                        end if;
                    end if;
                    if (hit_counter = 2) then 
                        ballcnt <= b"000";
                        hit_counter <= to_unsigned(0,3);
                        if (speed > 0) then 
                            speed <= speed -1;
                        end if;
                    end if;    
                    game <= moveBall;
                    rowcnt <= to_unsigned(1,6); --1,6
               when moveBall =>
                if ((ball_col = 0) or (ball_col = 79)) then 
                   ball_col <= to_unsigned(40, 7);
                   ball_row <= to_unsigned(30,6);
                elsif (ballcnt < speed) then 
                    ballcnt <= ballcnt + 1;
                else
                    ballcnt <= b"000";
                    if (y_count < y_thresh) then 
                        if ((y_change < 0) and (ball_row >1)) then 
                            ball_row <= ball_row - 1;
                        elsif ((y_change > 0) and (ball_row < 58)) then 
                            ball_row <= ball_row + 1;
                        elsif (ball_row = 1) then
                            ball_row <= ball_row + 1;
                            y_change <= 1;
                        elsif (ball_row = 58) then 
                            ball_row <= ball_row -1;
                            y_change <= -1;
                        end if;
                        y_count <= y_count + 1;                
                    else 
                        if (ball_col < 79) and (ball_col > 0) then 
                            ball_col <= ball_col + ball_dx;
                        end if;
                        y_count <= to_unsigned(0,4);
                    end if;
               end if;
               game <= addBall;
               when addPaddle =>
                    if (rowcnt < 59) then    --59
                            if ((rowcnt >= pad1_top) and (rowcnt < pad1_bot)) then
                                pong_bg(to_integer(rowcnt), 0) <= b"11";
                            else
                                if (ball_col /= 0) then 
                                    pong_bg(to_integer(rowcnt), 0) <= b"00";
                                else 
                                    if (ball_row /= rowcnt) then
                                        pong_bg(to_integer(rowcnt), 0) <= b"00";
                                    end if;
                                end if;
                            end if;
                            if ((rowcnt >= pad2_top) and (rowcnt < pad2_bot)) then 
                                pong_bg(to_integer(rowcnt), 79) <= b"11"; --79
                            else 
                                if (ball_col /= 79) then  --79
                                    pong_bg(to_integer(rowcnt), 79) <= b"00"; --79
                                else
                                    if (ball_row /= rowcnt) then
                                        pong_bg(to_integer(rowcnt), 79) <= b"00";--79
                                    end if;
                                end if;
                            end if;
                            rowcnt <= rowcnt + 1;
                   else 
                        if ((score1 = '1') or (score2 = '1')) then 
                            game <= startGame;
                            if (score1 = '1') then 
                                ball_col <= to_unsigned(76, 7);--76,7
                            else 
                                ball_col <= to_unsigned(3, 7); --3, 7
                            end if;      
                        else 
                            game <= idle;
                        end if;
                        rowcnt <= to_unsigned(0,6); --6
                        colcnt <= to_unsigned(0,7);--7
                    end if;
               when addBall =>
                    if (rowcnt < 59) then  --59
                        if (colcnt < 80) then --80
                            if (ball_row = rowcnt) and (ball_col = colcnt) then
                                pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"11";
                            else 
                                if (colcnt /= 39) and (colcnt /=40)  then --39, 40
                                  pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= b"00";
                                end if;
                           end if;
                           colcnt <= colcnt + 1;
                        else 
                            colcnt <= to_unsigned(0,7);--7
                            rowcnt <= rowcnt + 1;
                        end if;
                   else 
                        rowcnt <= to_unsigned(1,6); --6
                        game <= addPaddle;
                   end if;
               when others =>
                    null;
              end case;
             end if;
             end process;
             
	
	

	       
	
	       

end arch;