$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))

if ENV['COVERAGE'] 
  require 'simplecov'
  SimpleCov.start do
    add_filter '/spec/'
    add_filter '/.bundle/'
  end
end

require 'rspec'
require 'align'


# Requires supporting files with custom matchers and macros, etc,
# in ./support/ and its subdirectories.
Dir["#{File.dirname(__FILE__)}/support/**/*.rb"].each {|f| require f}

RSpec.configure do |config|
  
end
