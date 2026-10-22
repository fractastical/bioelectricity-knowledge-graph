"""
Based on code from Amitav Nott.
Two-channel, sweep-style oscilloscope.
"""

import sys, argparse, time, serial, math
import matplotlib.pyplot as plt
import numpy as np

# ==== CONFIGURATION ====
DEFAULT_BAUDRATE = 400000	# For communication from Nucleo board => host

# PLOT_WINDOW_SEC is the width of the X axis in seconds.
PLOT_WINDOW_SEC = 5

# READS_PER_SEC is our expectation of how often the Nucleo board will send us
# samples. We mostly care because this drives buffer sizes for our plots.
# More subtly, after each loop we sleep just long enough to wake up at about
# the next data burst.
READS_PER_SEC = 500

BUFFER_SIZE = math.ceil(READS_PER_SEC * PLOT_WINDOW_SEC)
RECORDING_SIZE = READS_PER_SEC * 10	# Dump at most 10 seconds of data

# Pause/unpause when the user hits 'p'.
g_paused = False

# Record the incoming data and save it to a file. The strategy:
# - Class Recording is a static class; we never bother to allocate any objects.
# - It allocates a two-channel data buffer recording_array[RECORDING_SIZE,2]
# - Assume that somebody calls Recording.save_pair() for every incoming data
#   point. Typically, save_pair() does nothing. But when we're in recording
#   mode (because the user typed an 'r'), we save the data into recording_array.
# - When recording_array[] fills up, we call Recording.dump() to dump it. We
#   always just dump to scope_out.csv; if they hit 'r' twice, the second time
#   overwrites the first.
class Recording:
    recording_array = np.zeros ([RECORDING_SIZE,2], dtype=np.int16)
    record_idx = 0
    is_recording = False

    # Save_pair() is a noop until somebody hits "r"; then it saves data into our
    # buffer until the buffer fills, and then says "Done".
    @classmethod
    def save_pair (cls, ch1, ch2):
        if (cls.is_recording):
            cls.recording_array[cls.record_idx,0] = ch1 * 4095
            cls.recording_array[cls.record_idx,1] = ch2 * 4095
            cls.record_idx += 1
        if (cls.record_idx == RECORDING_SIZE) and cls.is_recording:
            cls.is_recording = False
            print ("Done recording")
            cls.dump("scope_out.csv")
    @classmethod
    def dump (cls,filename):
        print (f"Dumping to {filename}")
        if (cls.record_idx == 0):
            print ("No saved data to dump to",filename)
            return
        with open (filename, "w") as outfile:
            for i in range (cls.record_idx):
                outfile.write (f"{cls.recording_array[i,0]},"
                               f"{cls.recording_array[i,1]}\n")
        print ("File dump is complete")
        cls.record_idx = 0
    @classmethod
    def start_recording(cls):
        print ("Recording...")
        cls.is_recording = True

# Parse the command line and return a data structure with the user input.
def get_command_args():
    parser = argparse.ArgumentParser(description='Python sweep-style scope')
    parser.add_argument('--stdin', action='store_true',
                        help='Read from stdin instead of serial')
    parser.add_argument('--port', type=str, help='Serial port name')
    parser.add_argument('--baud', type=int, default=DEFAULT_BAUDRATE,
                        help='Serial baudrate')
    return parser.parse_args()

# Figure out which serial port to read from, and return it as an ascii string
# (e.g., "COM3"). The hard work is done by serial.tools.list_ports.comports(),
# which returns a list of all possible ports.
def get_serial_port(args):
    import serial.tools.list_ports
    ports = serial.tools.list_ports.comports()
    port_names = [p.device for p in ports]
    if (args.port):		# User has specified a port
        if (not (args.port.upper() in port_names)):
         print(f"{args.port} isn't one of {port_names}, even case-insensitive!")
        return (args.port)
    if (len(ports)==1):	# Only one port is possible
        return (ports[0].device)
    if (len(port_names) == 0):
        print("You have no serial port connected; is your USB cord plugged in?")
        quit()
    assert len(port_names) > 1
    print (f"You have multiple serial input ports {port_names}. Please choose:")
    print ("	0. Exit")
    for i in range (len (port_names)):
        print (f"\t{i+1}: {port_names[i]}")
    numb = int (input("Number: "))
    if (numb==0):
        quit()
    return (port_names[numb-1])

# We use Matplotlib's mpl_connect() to call this handler for any key presses
# when our plot has the user's focus.
def keypress_handler(event):
    global g_paused
    # print (f"In handler for key {event.key}")

    if (event.key=='r'):
        Recording.start_recording()
    if (event.key=='p'):
        g_paused = not g_paused

def main():
    args = get_command_args()
    # Choose input source
    if args.stdin:
        input_source = sys.stdin.buffer
        print("reading input from stdin")
    else:
        port = get_serial_port(args)
        print (f"Using serial port {port}")
        ser = serial.Serial(port, args.baud, timeout=0)
        input_source = ser

    # Set up plotting in non-interactive mode. This means that plotting commands
    # don't take effect until we force redraw with .draw_idle() later.
    plt.ioff()
    fig, ax = plt.subplots()

    # Build a 1D array of the evenly-spaced X-axis points. We'll modify it each
    # time the scope trace sweeps to the left. E.g., it will start as [0...5],
    # then become [5...10], then [10...15], etc.
    x = np.linspace(0, PLOT_WINDOW_SEC, BUFFER_SIZE)

    # Build the two Y-value buffers and turn them into lines.
    # They're floats in [0,1] (initially all zeroes).
    buffer1 = np.zeros(BUFFER_SIZE, dtype=np.float32)	# Channel 1 = Lead 1
    buffer2 = np.zeros(BUFFER_SIZE, dtype=np.float32)	# Channel 2 = Lead 2
    line1, = ax.plot(x, buffer1, label="Lead 1")
    line2, = ax.plot(x, buffer2, label="Lead 2")

    # The cursor is a moving full-plot-height vertical red line . Incoming new
    # data gets drawn at the cursor.
    cursor_line = ax.axvline(x[0], color="red", linestyle="--")
    cursor = 0	# Counts circularly from [0...BUFFER_SIZE).

    # A bit more axis/label setup.
    ax.set_xlim(0, PLOT_WINDOW_SEC)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('ECG signal')
    ax.set_title('Python two-channel oscilloscope')
    ax.legend(loc="upper right")

    # For computing frames/sec and data_points-pairs/sec to add to the title.
    last_time = time.time()
    loops_this_sec = 0
    points_this_sec = 0

    # Update interval is our expected time between getting new points (in sec).
    update_interval = 1.0 / READS_PER_SEC
    leftover = b''

    # Set up recording/file dumping.
    handler=fig.canvas.mpl_connect('key_press_event', keypress_handler)

    print ("Type 'r' to start recording (data saved to scope_out.csv) or 'p' to pause/unpause the plot")
    print ("Exit by 'x'-ing out the plot window (not by ctrl-c)")
    try:
        while plt.fignum_exists(fig.number):	# Until the user 'x'es it out.
            loop_start = time.time()	# To sleep until the next expected data

            # Read any raw bytes
            if args.stdin:
                if (hasattr(input_source, "read1")):
                    raw = input_source.read1(1024)
                else:
                    raw = input_source.read(1024)
            else:
                raw = input_source.read(input_source.in_waiting or 1)

            if raw:
                # print (f"Read {len(raw)} bytes")
                # Prepend any leftover from previous read
                raw = leftover + raw

                # We expect 3-byte sequences with MSB markers:
                # [1|data2][0|data1][0|data0] -> 21 bits of payload (7 bits per
                # byte). Assemble the 21-bit payload and split into an 11-bit
                # channel and a 10-bit channel.
                i = 0			# Current index into raw[].
                last_processed = 0
                raw_len = len(raw)
                while i <= raw_len - 3:
                    b0 = raw[i]
                    b1 = raw[i+1]
                    b2 = raw[i+2]
                    # Start marker: first byte MSB=1, next two bytes MSB=0
                    if (b0 & 0x80) and not (b1 & 0x80) and not (b2 & 0x80):
                        b0 &= 0x7F
                        b1 &= 0x7F
                        b2 &= 0x7F
                        combined = (b0 << 14) | (b1 << 7) | b2  # 21-bit value
                        ch1 = (combined >> 10) & 0x7FF  # 11 MSB
                        ch2 = combined & 0x3FF          # 10 LSB
                        # Scale channels to 0-255 for plotting. // is division
                        # with a "floor" afterwards.
                        v1 = (ch1 * 255) / 0x7FF # since it was 11 bits
                        v2 = (ch2 * 255) / 0x3FF # since it was 10 bits
                        v1 = ch1 / 0x7FF # since it was 11 bits
                        v2 = ch2 / 0x3FF # since it was 10 bits
                        # Update the 1D Y-value arrays for plotting.
                        if (not g_paused):
                            buffer1[cursor] = v1
                            buffer2[cursor] = v2
                            cursor += 1
                        Recording.save_pair (v1, v2)
                        #print (f"Bytes {b0:03x}, {b1:03x}, {b2:03x}" \
                        #     + f"-> {combined:06x} -> {ch1:03x}, {ch2:03x}"\
                        #     + f"-> {v1}, {v2}")

                        # Wrap cursor if needed; change X axis when we wrap.
                        if cursor >= BUFFER_SIZE:
                            cursor = 0
                            x += PLOT_WINDOW_SEC
                            ax.set_xlim(x[0], x[-1])
                            line1.set_xdata(x)
                            line2.set_xdata(x)
                        # Suck away the bytes we just processed.
                        i += 3
                        last_processed = i
                        points_this_sec += 1	# To print points/sec.
                    else:
                        # Not a valid aligned triple -- skip the current byte
                        # and look for a start marker at the next byte. This
                        # typically happens at startup, and if the USB link ever
                        # loses a byte.
                        i += 1

                # We've searched through this entire set of new bytes. Keep any
                # remaining bytes (incomplete sequence) for the next read.
                leftover = raw[last_processed:]
                # Prevent leftover growing unbounded; at most 2 bytes are
                # needed to complete a 3-byte sequence
                if len(leftover) > 2:
                    leftover = leftover[-2:]
                    # print (f"Discarded {len(leftover)-2}B")

            # Redraw the updated traces.
            line1.set_ydata(buffer1)
            line2.set_ydata(buffer2)

            # Redraw the cursor in its new location.
            cursor_time = x[cursor]
            cursor_line.set_xdata([cursor_time, cursor_time])

            # Matplotlib magic !!???
            fig.canvas.draw_idle()	# Redraw when event loop is idle
            plt.pause(0.001)		# Allows the event loop to process.

            # Compute and display the frame rate. Essentially, the number of
            # read/process/display_update loops/second. We even count loops that
            # don't read anything!
            loops_this_sec += 1
            now = time.time()
            if now - last_time > 1.0:	# More than 1 second...
                fps = loops_this_sec / (now - last_time)
                pps = points_this_sec/ (now - last_time)
                ax.set_title(f'Python scope: {fps:.1f} FPS, {pps:.1f}pt/sec')
                last_time = now
                loops_this_sec = 0
                points_this_sec = 0

            # Update interval is our expected time between getting new points
            # (in sec). We'll do a sleep so as to roughly only run the outer
            # loop when we think data is ready. I don't think this is useful
            # except for niceness.
            elapsed = time.time() - loop_start
            if elapsed < update_interval:
                time.sleep(update_interval - elapsed)

    finally:
        print ("All done")
        if not args.stdin:
            input_source.close()
        fig.canvas.mpl_disconnect (handler)

if __name__ == "__main__":
    main()
