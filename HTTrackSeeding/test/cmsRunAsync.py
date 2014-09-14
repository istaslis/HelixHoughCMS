import sys
import os.path

print os.getcwd()
sys.path.append(os.getcwd())

for i, modulepath in enumerate(sys.argv):
    if i==0 : continue
    print 'module! '+str(i)+' ' + modulepath

    modulename = os.path.splitext(modulepath)[0]
    
    configfile = open(modulepath).read()
    
    module = __import__(modulename, fromlist=['outputfilename'])
    outname = os.path.splitext(module.outputfilename)[0]
    print module.outputfilename+' '+outname
    output1 = outname+'output1'
    output2 = outname+'output2'
    
    command = 'cmsRun '+modulepath+' 1>'+output1+' 2>'+output2
    print command
    
    
    exitcode = os.system(command)
        
#send an email!
    import smtplib

    # Import the email modules we'll need
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
    from email.MIMEBase import MIMEBase
    from email import Encoders
    
    msg = MIMEMultipart()
    
    msg.attach( MIMEText('Job is finished with exit code '+str(exitcode)+'! \n Command: '+command+'\n'+'Output file: '+module.outputfilename+'\n'+configfile+'\n') )
    files = [output1,output2]
    
    for f in files:
        if os.path.isfile(f) and os.path.getsize(f)<1000000:
            part = MIMEBase('application', "octet-stream")
            part.set_payload( open(f,"rb").read() )
            Encoders.encode_base64(part)
            part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(f))
            msg.attach(part)
            
            
            
    sender = os.environ['USER']+'@'+os.environ['HOSTNAME']
    receiver = 'lisniak@llr.in2p3.fr'
            
    # me == the sender's email address
    # you == the recipient's email address
    msg['Subject'] = 'Job is done!'
    msg['From'] = sender
    msg['To'] = receiver
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(sender, receiver, msg.as_string())
    s.quit()
            
