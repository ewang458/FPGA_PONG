library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;
library UNISIM;
use UNISIM.vcomponents.all;

entity pong is
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
end pong;

architecture arch of pong is

component debounce 
    port(
        clk: in std_logic;
        btn: in std_logic;
        btn_out: out std_logic
        );
end component;  
component score
    port(
        clk: in std_logic;
        score_flag  : in  std_logic;
        identity : in std_logic;
        hcount:   in unsigned(9 downto 0);
        vcount:  in unsigned(9 downto 0);
        char_pixel : out std_logic;
        current_score : out integer range 0 to 20;
        rst: in std_logic
        );
end component; 
	signal clkfb:    std_logic;
	signal clkfx:    std_logic;
	signal hcount:   unsigned(9 downto 0);
	signal vcount:   unsigned(9 downto 0);
	signal blank:    std_logic;
	signal frame:    std_logic;
	signal frame_prev: std_logic;
	signal pad2_top: unsigned(5 downto 0) := to_unsigned(23, 6);
	signal pad2_bot: unsigned(5 downto 0) := to_unsigned(35, 6);
	signal pad1_top: unsigned(5 downto 0) := to_unsigned(23, 6);
	signal pad1_bot: unsigned(5 downto 0) := to_unsigned(35, 6);
	signal ball_dx: signed(3 downto 0) := to_signed(2,4);
	signal ball_dy: signed(3 downto 0) := to_signed(2,4);
	signal y_change: signed(3 downto 0); 
	signal y_thresh: signed(3 downto 0);
	signal y_count: unsigned(3 downto 0);
	signal ball: std_logic;
	signal ball_y_top: unsigned(9 downto 0):= to_unsigned(240, 10);
	signal ball_y_bot: unsigned(9 downto 0) := to_unsigned(248,10);
	signal ball_x_right: unsigned (9 downto 0) := to_unsigned(328,10);
	signal ball_x_left: unsigned (9 downto 0) := to_unsigned(320,10);
	signal ballcnt: unsigned(2 downto 0);
	signal score1: std_logic := '0'; --left player
	signal score2: std_logic := '0'; --right player
	signal hit_counter: unsigned(2 downto 0);
	signal speed: unsigned(2 downto 0) := to_unsigned(4,3);
	signal up_1: std_logic;
	signal up_2: std_logic;
	signal down_1: std_logic;
	signal down_2: std_logic;
	signal rst: std_logic;
	type bg is array(59 downto 0,79  downto 0) of std_logic;
	signal pong_bg: bg;
	signal rowcnt: unsigned(5 downto 0) := (others => '0');
	signal colcnt : unsigned(6 downto 0) := (others => '0');
	signal pixel_jump: integer := 2;

	--signals for scorekeeping:
	signal score_left: integer range 0 to 11 := 0;
    signal char_pixel_left: std_logic;
    signal score_right: integer range 0 to 11 := 0;
    signal char_pixel_right: std_logic;
    signal score_flag: std_logic := '0';
    signal start_flag: std_logic := '0';
    signal left_id: std_logic := '0'; --left player identified by score.vhd with 0
    signal right_id: std_logic := '1'; --right player identified by score.vhd with 1
    signal reset_game: std_logic := '0'; --game reset caused by score reaching 11 instead of btn1 being pressed
	type FSM is (startFrame, startGame, idle, initialize, updatePaddle, updateBall, moveBall, addPaddle, addBall);
	signal game: FSM;

begin
    up1debounce: debounce port map(clk => clkfx, btn => up1, btn_out => up_1);
    up2debounce: debounce port map(clk => clkfx, btn => up2, btn_out => up_2);
    down1debounce: debounce port map(clk => clkfx, btn => down1, btn_out => down_1);
    down2debounce: debounce port map(clk => clkfx, btn => down2, btn_out => down_2);
    left_player: score port map(
                clk => clkfx, 
                score_flag => score1, 
                identity => left_id, 
                hcount => hcount, 
                vcount => vcount, 
                char_pixel => char_pixel_left, 
                current_score => score_left,
                rst => reset
            );
    right_player: score port map(
                clk => clkfx, 
                score_flag => score2, 
                identity => right_id, 
                hcount => hcount, 
                vcount => vcount, 
                char_pixel => char_pixel_right, 
                current_score => score_right,
                rst => reset
            );
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
	
	red<=b"00" when blank='1' else b"11" when ball = '1' else
       b"11" when (char_pixel_left = '1' or char_pixel_right = '1') else
       b"11" when (((vcount(9 downto 3) < 60) and (hcount(9 downto 3) < 80))and (pong_bg(to_integer(vcount(9 downto 3)), to_integer(hcount(9 downto 3))) = '1') )else
       b"00";
	green<=b"00" when blank='1' else b"00" when ball = '1' else
       b"11" when (char_pixel_left = '1' or char_pixel_right = '1') else
       b"11" when (((vcount(9 downto 3) < 60) and (hcount(9 downto 3) < 80))and (pong_bg(to_integer(vcount(9 downto 3)), to_integer(hcount(9 downto 3))) = '1') )else
       b"00";
	blue<=b"00" when blank='1' else b"00" when ball = '1' else
       b"11" when (char_pixel_left = '1' or char_pixel_right = '1') else
       b"11" when (((vcount(9 downto 3) < 60) and (hcount(9 downto 3) < 80))and (pong_bg(to_integer(vcount(9 downto 3)), to_integer(hcount(9 downto 3))) = '1') )else
       b"00";
	--OUTPUT TO VGA (SCALE DOWN PONG_BG ARRAY BY DIVIDING BY 8, and ensure hcount within blanking range);\
	
	process(hcount, vcount, reset, start_flag)
	begin
	if (reset = '1') then 
	    ball <= '0';	   
	elsif (start_flag = '1') then  
	   if ((hcount >= ball_x_left) and (hcount < ball_x_right)) and
	       ((vcount >= ball_y_top) and (vcount < ball_y_bot))then     
	           ball<= '1';
	   else 
               ball<= '0';
	   end if;
	else 
	   ball <= '0';
	end if;
	end process;
	
	process(clkfx, reset)
	begin 
	if (reset = '1') then 
	   pixel_jump <= 2;
	   start_flag <= '0';
	   ballcnt <= b"000";
	   ball_dx <= to_signed(2,4);
       ball_dy <= to_signed(2,4);
       speed <= to_unsigned(4,3);
       ball_y_top <= to_unsigned(240, 10);
	   ball_y_bot <= to_unsigned(248,10);
	   ball_x_right <= to_unsigned(328,10);
	   ball_x_left <= to_unsigned(320,10);
       score1 <= '0'; --reset score flags here for now, may have to be dealt with differently when score added
       score2 <= '0'; 
       ballcnt <= to_unsigned(0,3);
       hit_counter <= to_unsigned(0,3);
	elsif rising_edge(clkfx) then
	frame_prev <= frame;
	if ((start = '1') and (game = startGame)) then 
	   pixel_jump <= 2;
	   hit_counter <= to_unsigned(0,3);
       speed <= to_unsigned(4,3);
       score1 <= '0'; --reset score flags here for now, may have to be dealt with differently when score added
       score2 <= '0';  
       start_flag <= '1';
    end if;
	if (frame = '1') and (frame_prev = '0') then 
	if(start_flag = '1') then
	   if (ballcnt < 0) then  --counter for changing ball speed
            ballcnt <= ballcnt + 1;
       else 
            ballcnt <= b"000"; --reset counter
            if (abs(ball_dy) = pixel_jump) then -- if  y ball speed is 1 then we do diagonal movement 
                y_count <= b"0000"; --reset y counter for L movement 
                if ((ball_x_left = 8) or (ball_x_right = 632)) then --if ball in paddle range
                    if (ball_x_left = 8) then  -- if on left side 
                        if ((ball_y_top >= (pad1_top&"000")) and (ball_y_top < ((pad1_top + 3)&"000"))) or ((ball_y_bot >= (pad1_top&"000")) and (ball_y_bot < ((pad1_top + 3)&"000"))) then  -- if on top 3 pixels on paddle hit
                            ball_x_left <= ball_x_left + pixel_jump; --change x valueand velcoity
                            ball_x_right <= ball_x_right + pixel_jump;
                            ball_dx <= to_signed(pixel_jump,4); 
                        if (ball_dy < 0) then  --if ball was moving up, increase angle and keep moving up
                            ball_dy <= ball_dy -pixel_jump;
                            y_thresh <= abs(ball_dy - pixel_jump);
                            y_change <= to_signed(-pixel_jump,4);
                        else   --if ball was moving down, increase angle and now move up 
                            ball_dy <= -ball_dy -pixel_jump;
                            y_thresh <= abs(-ball_dy -pixel_jump);
                            y_change <= to_signed(-pixel_jump,4);
                        end if;
                        hit_counter <= hit_counter + 1; -- keep track of paddle hits
                    elsif ((ball_y_top >= ((pad1_top +3)&"000")) and (ball_y_top <((pad1_top +9)&"000"))) or ((ball_y_bot >= ((pad1_top +3)&"000")) and (ball_y_bot <((pad1_top +9)&"000")))then --if ball is in middle, decrease angle here
                        ball_x_left <= ball_x_left + pixel_jump; --change x valueand velcoity
                        ball_x_right <= ball_x_right + pixel_jump;
                        ball_dx <= to_signed(pixel_jump, 4);
                        ball_dy <= to_signed(0,4); --when y is 1 or -1 can only go to 0 horizontal movement 
                        y_thresh <= to_signed(0,4);
                        y_change <= to_signed(0,4);
                        hit_counter <= hit_counter + 1;
                    elsif ((ball_y_top >= ((pad1_top + 9)&"000")) and (ball_y_top < (pad1_bot&"000"))) or ((ball_y_bot >= ((pad1_top + 9)&"000")) and (ball_y_bot < (pad1_bot&"000"))) then --if bottom section of paddle
                        ball_x_left <= ball_x_left + pixel_jump; --change x valueand velcoity
                        ball_x_right <= ball_x_right + pixel_jump;
                        ball_dx <= to_signed(pixel_jump, 4);
                        if (ball_dy < 0) then --if ball was moving up, make it go down now increase angle
                            ball_dy <= -ball_dy + pixel_jump;
                            y_thresh <= -ball_dy + pixel_jump;
                            y_change <= to_signed(pixel_jump, 4);
                        else -- if ball was moving down, keep direction and increase angle 
                            ball_dy <= ball_dy + pixel_jump;
                            y_thresh <= ball_dy + pixel_jump;
                            y_change <= to_signed(pixel_jump,4);
                        end if;
                        hit_counter <= hit_counter + 1;
                    else --if not hitting paddle, continue in direction 
                        if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                                ball_x_right <= to_unsigned(639, 10);
                                ball_x_left <= to_unsigned(631, 10);
                         else
                             ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                             ball_x_right <= ball_x_right + to_integer(ball_dx);
                         end if;              
                         ball_y_top <= ball_y_top + to_integer(ball_dy); --change x valueand velcoity
                         ball_y_bot <= ball_y_bot + to_integer(ball_dy);
                    end if;
                elsif (ball_x_right = 632) then 
                    if ((ball_y_top >= (pad2_top&"000")) and (ball_y_top < ((pad2_top + 3)&"000"))) or ((ball_y_bot >= (pad2_top&"000")) and (ball_y_bot < ((pad2_top + 3)&"000"))) then  -- if on top 3 pixels on paddle hit
                        ball_x_left <= ball_x_left - pixel_jump; --change x valueand velcoity
                        ball_x_right <= ball_x_right - pixel_jump;
                        ball_dx <= to_signed(-pixel_jump,4); 
                        if (ball_dy < 0) then  --if ball was moving up, increase angle and keep moving up
                            ball_dy <= ball_dy -pixel_jump;
                            y_thresh <= abs(ball_dy - pixel_jump);
                            y_change <= to_signed(-pixel_jump,4);
                        else   --if ball was moving down, increase angle and now move up 
                            ball_dy <= -ball_dy -pixel_jump;
                            y_thresh <= abs(-ball_dy -pixel_jump);
                            y_change <= to_signed(-pixel_jump,4);
                        end if;
                        hit_counter <= hit_counter + 1; -- keep track of paddle hits
                    elsif ((ball_y_top >= ((pad2_top +3)&"000")) and (ball_y_top <((pad2_top +9)&"000"))) or ((ball_y_bot >= ((pad2_top +3)&"000")) and (ball_y_bot <((pad2_top +9)&"000"))) then --if ball is in middle, decrease angle here
                        ball_x_left <= ball_x_left - pixel_jump; --change x valueand velcoity
                        ball_x_right <= ball_x_right - pixel_jump;
                        ball_dx <= to_signed(-pixel_jump, 4);
                        ball_dy <= to_signed(0,4); --when y is 1 or -1 can only go to 0 horizontal movement 
                        y_thresh <= to_signed(0,4);
                        y_change <= to_signed(0,4);
                        hit_counter <= hit_counter + 1;
                    elsif ((ball_y_top >= ((pad2_top + 9)&"000")) and (ball_y_top < (pad2_bot&"000"))) or ((ball_y_bot >= ((pad2_top + 9)&"000")) and (ball_y_bot < (pad2_bot&"000"))) then --if bottom section of paddle
                        ball_x_left <= ball_x_left - pixel_jump; --change x valueand velcoity
                        ball_x_right <= ball_x_right - pixel_jump;
                        ball_dx <= to_signed(-pixel_jump, 4);
                        if (ball_dy < 0) then --if ball was moving up, make it go down now increase angle
                            ball_dy <= -ball_dy + pixel_jump;
                            y_thresh <= -ball_dy + pixel_jump;
                            y_change <= to_signed(2, 4);
                        else -- if ball was moving down, keep direction and increase angle 
                            ball_dy <= ball_dy + 2;
                            y_thresh <= ball_dy + 2;
                            y_change <= to_signed(2,4);
                        end if;
                        hit_counter <= hit_counter + 1;
                     else --if not hitting paddle, continue in direction 
                        if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                                ball_x_right <= to_unsigned(639, 10);
                                ball_x_left <= to_unsigned(631, 10);
                         else
                             ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                             ball_x_right <= ball_x_right + to_integer(ball_dx);
                         end if;              
                         ball_y_top <= ball_y_top + to_integer(ball_dy); --change x valueand velcoity
                         ball_y_bot <= ball_y_bot + to_integer(ball_dy);
                     end if;
                 end if; 
              elsif (ball_x_left = 0) then --if ball hits left side of screen (player 2 score)
                        ball_dx <= to_signed(2,4); --set x velcoity so that it si a serve from player 1
                        score2 <= '1'; --player 2 (right) score
                        ball_x_right <= to_unsigned(32,10);
	                    ball_x_left <= to_unsigned(24,10);
                        start_flag <= '0';
                  elsif (ball_x_right = 639) then  --if ball hits right side of screen ( player 1 score)
                        ball_dx <= to_signed(-2,4); --plater 2 serves 
                        score1 <= '1';   --player 1 (left) score  
                        ball_x_right <= to_unsigned(616,10);
	                    ball_x_left <= to_unsigned(608,10);    
                        start_flag <= '0';               
                  elsif (ball_y_bot = 472) then --if hits bottom of screen, bounce 
                        ball_y_bot <= ball_y_bot -2;
                        ball_y_top <= ball_y_top -2;
                        ball_dy <= to_signed(-2, 4);
                        if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                                ball_x_right <= to_unsigned(639, 10);
                                ball_x_left <= to_unsigned(631, 10);
                         else
                             ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                             ball_x_right <= ball_x_right + to_integer(ball_dx);
                         end if;
                  elsif (ball_y_top = 8) then  --if hits top of screen bounce 
                        ball_y_bot <= ball_y_bot + 2;
                        ball_y_top <= ball_y_top + 2;
                        ball_dy <= to_signed(2, 4);
                        if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                             ball_x_right <= to_unsigned(639, 10);
                             ball_x_left <= to_unsigned(631, 10);
                        else
                             ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                             ball_x_right <= ball_x_right + to_integer(ball_dx);
                        end if;
                  else  -- otherwise move diagonally 
                       if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                             ball_x_right <= to_unsigned(639, 10);
                             ball_x_left <= to_unsigned(631, 10);
                        else
                             ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                             ball_x_right <= ball_x_right + to_integer(ball_dx);
                        end if;
                        ball_y_top <= ball_y_top + to_integer(ball_dy); --change x valueand velcoity
                        ball_y_bot <= ball_y_bot + to_integer(ball_dy);
                  end if; 
           else  -- if y velcoity is not equal to 1
               if ((ball_x_left = 8) or (ball_x_right = 632)) then --hitting paddles 
                    if (ball_x_left = 8) then -- if hits player 1 paddle 
                        if ((ball_y_top >= (pad1_top&"000")) and (ball_y_top < ((pad1_top + 3)&"000"))) or ((ball_y_bot >= (pad1_top&"000")) and (ball_y_bot < ((pad1_top + 3)&"000"))) then  -- hits top section of paddle
                                y_count <= b"0000"; --set counter down to 0
                                ball_x_left <= ball_x_left + 2;
                                ball_x_right <= ball_x_right + 2; -- move away from paddle so that we don't have weird logic 
                                ball_dx <= to_signed(2,4);
                                if (ball_dy < 0) then --if ball is  moving up then keep moving up and then increase angle
                                    if ((ball_dy - 2) < -8) then --cap angle at y =4
                                        ball_dy <= to_signed(-8,4);
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(-2,4);
                                    else  
                                        ball_dy <= ball_dy -2;
                                        y_thresh <= abs(ball_dy - 2);
                                        y_change <= to_signed(-2,4);
                                        
                                    end if;
                                else  --if was going down make it go up and increase angle 
                                    if ((-ball_dy - 2) < -8) then --cap angle at y = 4 
                                        ball_dy <= to_signed(-8,4);
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(-2,4);
                                        
                                    else --if was going up keep up increaes angle 
                                        ball_dy <= -ball_dy -2;
                                        y_thresh <= abs(-ball_dy -2);
                                        y_change <= to_signed(-2,4);
                                   end if;
                                end if;
                                hit_counter <= hit_counter + 1;
                           elsif ((ball_y_top >= ((pad1_top +3)&"000")) and (ball_y_top <((pad1_top +9)&"000"))) or ((ball_y_bot >= ((pad1_top +3)&"000")) and (ball_y_bot <((pad1_top +9)&"000"))) then --if middle paddle
                                y_count <= b"0000"; --reset counter
                                ball_x_left <= ball_x_left + 2;
                                ball_x_right <= ball_x_right + 2; -- move away from paddle so that we don't have weird logic 
                                ball_dx <= to_signed(2, 4); 
                                if (ball_dy < 0) then --decrease angle
                                    ball_dy <= ball_dy  + 2;  
                                    y_thresh <= abs(ball_dy + 2);
                                    y_change <= to_signed(-2,4);
                                else --if angle was 0 then move downwards, otherwise decreasing angle for moving down
                                    ball_dy <= ball_dy -2;
                                    y_thresh <= ball_dy - 2;
                                    y_change <= to_signed(2,4);
                                end if;
                                hit_counter <= hit_counter + 1;
                           elsif ((ball_y_top >= ((pad1_top + 9)&"000")) and (ball_y_top < (pad1_bot&"000"))) or ((ball_y_bot >= ((pad1_top + 9)&"000")) and (ball_y_bot < (pad1_bot&"000"))) then --if on bottom of paddle 
                                y_count <= b"0000";
                                ball_x_left <= ball_x_left + 2;
                                ball_x_right <= ball_x_right + 2; -- move away from paddle so that we don't have weird logic 
                                ball_dx <= to_signed(2, 4);
                                if (ball_dy < 0) then -- if ball was moving up then increase angle keep moving up
                                    if (-ball_dy + 2 > 8) then 
                                        ball_dy <= to_signed(8,4); --cap y at 4
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(2,4);
                                    else  -- if ball was moving down then move up and increaes angle 
                                        ball_dy <= -ball_dy + 2; 
                                        y_thresh <= -ball_dy + 2;
                                        y_change <= to_signed(2, 4);
                                    end if;
                                elsif (ball_dy > 0) then --if ball was moving down then keep and increase angle
                                    if (ball_dy + 2 > 8) then 
                                        ball_dy <= to_signed(8,4);
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(2,4);
                                    else
                                        ball_dy <= ball_dy + 2;
                                        y_thresh <= ball_dy + 2;
                                        y_change <= to_signed(2,4);
                                    end if;
                                else -- if ball was zero move it down 
                                    ball_dy <= to_signed(2,4);
                                    y_thresh <= to_signed(2,4);
                                    y_change <= to_signed(2,4);
                               end if;
                                hit_counter <= hit_counter + 1;
                           else --if not hit paddle keep mvoing as planned
                                if (y_count < unsigned(y_thresh)) then 
                                    ball_y_top <= ball_y_top+ to_integer(y_change);
                                    ball_y_bot <= ball_y_bot+ to_integer(y_change);
                                    y_count <= y_count + 2;
                                else 
                                    y_count <= b"0000";
                                    if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                                        ball_x_right <= to_unsigned(639, 10);
                                        ball_x_left <= to_unsigned(631, 10);
                                    else
                                        ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                                        ball_x_right <= ball_x_right + to_integer(ball_dx);
                                    end if;
                                end if;        
                           end if;
                    elsif (ball_x_right = 632) then  -- ITS THE SAME LOGIC FOR PADDLE 2 JUST LOOK AT PADDLE 1 COMMENTS 
                        if ((ball_y_top >= (pad2_top&"000")) and (ball_y_top < ((pad2_top + 3)&"000"))) or ((ball_y_bot >= (pad2_top&"000")) and (ball_y_bot < ((pad2_top + 3)&"000"))) then  -- hits top section of paddle
                                y_count <= b"0000"; --set counter down to 0
                                ball_x_left <= ball_x_left - 2;
                                ball_x_right <= ball_x_right - 2; -- move away from paddle so that we don't have weird logic 
                                ball_dx <= to_signed(-2,4);
                                if (ball_dy < 0) then --if ball is  moving up then keep moving up and then increase angle
                                    if ((ball_dy - 2) < -8) then --cap angle at y =4
                                        ball_dy <= to_signed(-8,4);
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(-2,4);
                                    else  
                                        ball_dy <= ball_dy -2;
                                        y_thresh <= abs(ball_dy - 2);
                                        y_change <= to_signed(-2,4);
                                        
                                    end if;
                                else  --if was going down make it go up and increase angle 
                                    if ((-ball_dy - 2) < -8) then --cap angle at y = 4 
                                        ball_dy <= to_signed(-8,4);
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(-2,4);
                                        
                                    else --if was going up keep up increaes angle 
                                        ball_dy <= -ball_dy -2;
                                        y_thresh <= abs(-ball_dy -2);
                                        y_change <= to_signed(-2,4);
                                   end if;
                                end if;
                                hit_counter <= hit_counter + 1;
                           elsif ((ball_y_top >= ((pad2_top +3)&"000")) and (ball_y_top <((pad2_top +9)&"000"))) or ((ball_y_bot >= ((pad2_top +3)&"000")) and (ball_y_bot <((pad2_top +9)&"000")))then --if middle paddle
                                y_count <= b"0000"; --reset counter 
                                ball_x_left <= ball_x_left - 2;
                                ball_x_right <= ball_x_right - 2; -- move away from paddle so that we don't have weird logic 
                                ball_dx <= to_signed(-2, 4); 
                                if (ball_dy < 0) then --decrease angle
                                    ball_dy <= ball_dy  + 2;  
                                    y_thresh <= abs(ball_dy + 2);
                                    y_change <= to_signed(-2,4);
                                else --if angle was 0 then move downwards, otherwise decreasing angle for moving down
                                    ball_dy <= ball_dy -2;
                                    y_thresh <= ball_dy - 2;
                                    y_change <= to_signed(2,4);
                                end if;
                                hit_counter <= hit_counter + 1;
                           elsif ((ball_y_top >= ((pad2_top + 9)&"000")) and (ball_y_top < (pad2_bot&"000"))) or ((ball_y_bot >= ((pad2_top + 9)&"000")) and (ball_y_bot < (pad2_bot&"000"))) then --if on bottom of paddle 
                                y_count <= b"0000";
                                ball_x_left <= ball_x_left - 2;
                                ball_x_right <= ball_x_right - 2; -- move away from paddle so that we don't have weird logic 
                                ball_dx <= to_signed(-2, 4);
                                if (ball_dy < 0) then -- if ball was moving up then increase angle keep moving up
                                    if (-ball_dy + 2 > 8) then 
                                        ball_dy <= to_signed(8,4); --cap y at 4
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(2,4);
                                    else  -- if ball was moving down then move up and increaes angle 
                                        ball_dy <= -ball_dy + 2; 
                                        y_thresh <= -ball_dy + 2;
                                        y_change <= to_signed(2, 4);
                                    end if;
                                elsif (ball_dy > 0) then --if ball was moving down then keep and increase angle
                                    if (ball_dy + 2 > 8) then 
                                        ball_dy <= to_signed(8,4);
                                        y_thresh <= to_signed(8,4);
                                        y_change <= to_signed(2,4);
                                    else
                                        ball_dy <= ball_dy + 2;
                                        y_thresh <= ball_dy + 2;
                                        y_change <= to_signed(2,4);
                                    end if;
                                else -- if ball was zero move it down 
                                    ball_dy <= to_signed(2,4);
                                    y_thresh <= to_signed(2,4);
                                    y_change <= to_signed(2,4);
                               end if;
                                hit_counter <= hit_counter + 1;
                           else --if not hit paddle keep mvoing as planned
                                if (y_count < unsigned(y_thresh)) then 
                                    ball_y_top <= ball_y_top+ to_integer(y_change);
                                    ball_y_bot <= ball_y_bot+ to_integer(y_change);
                                    y_count <= y_count + 2;
                                else 
                                    y_count <= b"0000";
                                    if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                                         ball_x_right <= to_unsigned(639, 10);
                                         ball_x_left <= to_unsigned(631, 10);
                                    else
                                         ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                                         ball_x_right <= ball_x_right + to_integer(ball_dx);
                                    end if;
                                end if;        
                           end if;
                       end if;
                  elsif (ball_x_left = 0) then  -- if ball is at edges increment score and say person scored 
                        ball_dx <= to_signed(2,4);
                        ball_dy <= y_change;
                        y_count <= b"0000";
                        score2 <= '1';
                        ball_x_right <= to_unsigned(32,10);
	                    ball_x_left <= to_unsigned(24,10);
                        start_flag <= '0';
                  elsif (ball_x_right = 639) then 
                        ball_dx <= to_signed(-2,4);
                        ball_dy <= y_change;
                        y_count <= b"0000";
                        score1 <= '1';
                        ball_x_right <= to_unsigned(616,10);
	                    ball_x_left <= to_unsigned(608,10);
                        start_flag <= '0';
                  elsif (ball_y_bot = 472) then  -- bounce if at top and bottom 
                        if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                             ball_x_right <= to_unsigned(639, 10);
                             ball_x_left <= to_unsigned(631, 10);
                        else
                             ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                             ball_x_right <= ball_x_right + to_integer(ball_dx);
                        end if;
                        ball_y_top <= ball_y_top - 2;
                        ball_y_bot <= ball_y_bot - 2;
                        y_change <= -y_change;
                        ball_dy <= -ball_dy;
                        y_count <= b"0000";
                  elsif (ball_y_top = 8) then 
                       if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                             ball_x_right <= to_unsigned(639, 10);
                             ball_x_left <= to_unsigned(631, 10);
                        else
                             ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                             ball_x_right <= ball_x_right + to_integer(ball_dx);
                        end if;
                        ball_y_top <= ball_y_top + 2;
                        ball_y_bot <= ball_y_bot + 2;
                        ball_dy <= -ball_dy;
                        y_change <= -y_change;
                        y_count <= b"0000";
                  else  --otherwise move as normal 
                        if (y_count < unsigned(y_thresh)) then 
                            ball_y_top <= ball_y_top+ to_integer(y_change);
                            ball_y_bot <= ball_y_bot+ to_integer(y_change);
                            y_count <= y_count + 2;
                        else 
                            y_count <= b"0000";
                            if ((ball_x_right + to_integer(ball_dx)) > 639) then  
                                 ball_x_right <= to_unsigned(639, 10);
                                 ball_x_left <= to_unsigned(631, 10);
                            else
                                 ball_x_left <= ball_x_left + to_integer(ball_dx); --change x valueand velcoity
                                 ball_x_right <= ball_x_right + to_integer(ball_dx);
                            end if;
                        end if;      
                  end if; 
                end if;                                      
                end if;    
                if (hit_counter = 2) then --if 2 hits then increase speed
                    ballcnt <= b"000";
                    hit_counter <= to_unsigned(0,3);
                    if (speed > 1) then 
                        speed <= speed -1;
                    end if;
                end if;  
                end if;
    end if;
    end if;
	end process;

	
    
    process(clkfx)
    begin
    score_flag <= '0';
    if (reset = '1' or reset_game = '1') then 
        --reset ball position and paddles (SCORE TOO) 
        pad1_top <= to_unsigned(23,6);
        pad1_bot <= to_unsigned(35,6);
        pad2_top <= to_unsigned(23,6);
        pad2_bot <= to_unsigned(35,6);
        rowcnt <= to_unsigned(0,6);
        colcnt <= to_unsigned(0,7);
        reset_game <= '0'; --a 1 here will prompt score to reset internally on both players in score module
        game <= startFrame; -- go to smaller frame 
        
    elsif (rising_edge(clkfx)) then
        if (reset_game = '1' and game = startFrame) then
            reset_game <= '0';
        elsif (score_left = 11 or score_right = 11) then
            reset_game <= '1';
        end if;
        case (game) is --add game start and game end states and reset capabilty 
            when startFrame =>  --draw game screen before start button state
                if (rowcnt < 60) then 
                    if (colcnt < 80) then 
                        if ((colcnt = 0) or (colcnt = 79)) then 
                            if ((rowcnt >= pad1_top) and (rowcnt < pad1_bot)) then 
                                pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= '1';
                            else
                                if ((rowcnt > 0) and (rowcnt < 59)) then
                                    pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= '0';
                                elsif ((rowcnt = 0) or (rowcnt = 59)) then 
                                    pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= '1';
                                end if;
                            end if;
                        elsif ((colcnt = 39) or (colcnt = 40)) then 
                            pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= '1';
                        elsif ((rowcnt = 0) or (rowcnt = 59)) then 
                            pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= '1';
                        else                     
                            pong_bg(to_integer(rowcnt), to_integer(colcnt)) <= '0';
                        end if;
                        colcnt <= colcnt + 1;
                    else 
                        colcnt <= to_unsigned(0,7);
                        rowcnt <= rowcnt + 1;
                    end if;
                else 
                    
                    game <= startGame;
                end if;                      
            when startGame => --wait for start button to be pushed, reset vlaues 
                if (start = '1') then 
                    game <= idle; --go to idle state and start game 
                end if;
            when idle => 
                if (frame = '1') then  --when vga enters blank period we start updating array
                    game <= updatePaddle;
                end if;
            when updatePaddle => --add paddle on screen, move with buttons
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
                rowcnt <= to_unsigned(1,6);
                game <= addPaddle;
               when addPaddle =>
                    if (rowcnt < 59) then    
                            if ((rowcnt >= pad1_top) and (rowcnt < pad1_bot)) then --draw paddle 1
                                pong_bg(to_integer(rowcnt), 0) <= '1';
                            else
                                pong_bg(to_integer(rowcnt), 0) <= '0';
                            end if; -- do the same thing 
                            if ((rowcnt >= pad2_top) and (rowcnt < pad2_bot)) then 
                                pong_bg(to_integer(rowcnt), 79) <= '1'; --79
                            else
                                pong_bg(to_integer(rowcnt), 79) <= '0'; --79
                            end if;
                            rowcnt <= rowcnt + 1;
                   else --after this go to startGame and wait for start button again if score occured
                        if ((score1 = '1') or (score2 = '1')) then
                            score_flag <= '1'; --score_flag -> '1' to prompt score change
                            --this is probably where we can include a flag for score to actually update!! 
                            game <= startGame;
                        else 
                            game <= idle; --otherwise keep moving throuhg, wiat for next frame
                        end if;
                        rowcnt <= to_unsigned(0,6); --6
                        colcnt <= to_unsigned(0,7);--7
                    end if;
               when others =>
                    null;
              end case;
             end if;
             end process;
end arch;
