import argparse
parser = argparse.ArgumentParser(description='test CHAP entrypoint script')
parser.add_argument('greeting', type = str, help = 'say hi')
args = parser.parse_args() 
if args.greeting == 'hi':
  print('yo whats up my man')
