tell application "Mail"

     set recipientName to "NaluCFD"
     set recipientAddress to "nalucfd@gmail.com"
     set theDate to date string of (current date)
     set theSubject to "Nalu Mac X Nightly Build and Test..." & theDate
     set buildNaluDocName to "Users:naluIt:gitHubWork:nightlyBuildAndTest:NaluBuild.txt"
     set buildTrilinosDocName to "Users:naluIt:gitHubWork:nightlyBuildAndTest:TrilinosBuild.txt"
     set rtestDocName to "Users:naluIt:gitHubWork:nightlyBuildAndTest:NaluRtest.txt"
     set theContent to "Nalu Build and Test Results Attached: " 
     theDate
 
        ##Create the message
        set theMessage to make new outgoing message with properties {subject:theSubject, content:theContent, visible:true}
 
        ##Set a recipient
        tell theMessage
                make new to recipient with properties {address:"naluCFD@gmail.com"}
                make new attachment with properties {file name:buildNaluDocName as alias}
                make new attachment with properties {file name:buildTrilinosDocName as alias}
                make new attachment with properties {file name:rtestDocName as alias}
                ##Send the Message
                send
        end tell
end tell
